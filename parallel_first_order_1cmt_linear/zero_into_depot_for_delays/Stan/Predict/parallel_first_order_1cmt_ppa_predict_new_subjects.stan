// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// This is flexible so that any number of parallel absorption processes can take 
//   place (n_depots >= 2).
// A delay in absorption for each process is implemented with a zero-order 
//   distributive delay (like an infusion into the depot) to bring about a
//   delayed absorption. The user can choose whether the fastest process has a 
//   delay by setting n_depots_with_delay to (n_depots - 1) for no delay or to
//   n_depots for a delay in the data input
// There will be n_depots KA values for each subject. i.e., TVKA and KAi are 
//   vectors. KA is intended to be ordered so that 
//   KA_1 < KA_2 < ... < KA_n_depots. This doesn't do this exactly, but it will 
//   have TVKA_1 < TVKA_2 < ... < TVKA_n_depots and hope the individual effects 
//   don't mess that up
// There will be n_depots_with_delay DUR values for each subject. i.e., TVDUR 
//   and DURi are vectors. DUR is intended to be ordered so that 
//   DUR_1 > DUR_2 > ... > DUR_n_depots_with_delay. If 
//   n_depots_with_delay = n_depots, then the fastest process has no delay. This
//   doesn't do this exactly, but it will have 
//   TVDUR_1 > TVDUR_2 > ... > TVDUR_(n_depots - 1) and hope the individual
//   effects don't mess that up
// TVFRAC is a simplex of length n_depots that tells how much of the dose goes 
//   into each absorption process. There is no IIV on FRAC
// IIV on CL, VC, KA, DUR (full covariance matrix)
// proportional plus additive error - DV = IPRED*(1 + eps_p) + eps_a
// General ODE solution using Torsten to get out individual estimates of AUC, 
//   Cmax, Tmax, ...
// Predictions are generated from a normal that is truncated below at 0

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  // An ODE that is written to take n_depots as an input. It should work for any
  //   n_depots >= 2
  vector parallel_1cmt_ode(real t, vector y, array[] real params, 
                           array[] real x_r, array[] int x_i){
    
    int n_depots = x_i[1];
    
    real cl = params[n_depots + 1];
    real vc = params[n_depots + 2];
    real ke = cl/vc;
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real slope = to_row_vector(params[1:n_depots])*y[1:n_depots] - 
                                                            ke*y[n_depots + 1];
    real x = slope > 0 && y[n_depots + 1]/vc > y[n_depots + 3] ? slope/vc : 0;
    real z = t <= t_1 || (slope > 0 && t >= t_1 && t <= t_2 && y[n_depots + 1]/vc > y[n_depots + 3]) ? 1 : 0;
    
    vector[n_depots + 4] dydt;

    for(i in 1:n_depots){
      dydt[i] = -params[i]*y[i];
    }
    dydt[n_depots + 1] = slope;
    dydt[n_depots + 2] = t >= t_1 && t <= t_2 ? y[n_depots + 1] : 0; // AUC_t_1-t_2
    dydt[n_depots + 3] = x;                                          // C_max
    dydt[n_depots + 4] = z;                                          // t_max
    
    return dydt;
    
  }
  
}
data{
  
  int n_subjects;
  int n_subjects_new;
  int n_time_new;
  array[n_time_new] real time;
  array[n_time_new] real amt;
  array[n_time_new] int cmt;
  array[n_time_new] int evid;
  array[n_time_new] real ii;
  array[n_time_new] int addl;
  array[n_time_new] int ss;
  array[n_subjects_new] int subj_start;
  array[n_subjects_new] int subj_end;
  
  int<lower = 2> n_depots; 
  int<lower = n_depots - 1, upper = n_depots> n_depots_with_delay;
  
  real<lower = 0> t_1;   // Time at which to start SS calculations (AUC_ss, C_max_ss, ...)
  real<lower = t_1> t_2; // Time at which to end SS calculations (AUC_ss, C_max_ss, ...)
 
}
transformed data{ 
  
  int n_random = n_depots + n_depots_with_delay + 2; // Number of random effects
  int n_cmt = n_depots + 4;                          // Number of states in the ODEs
  
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  array[1, 2] real x_r = {{t_1, t_2}};
  array[1, 1] int x_i = {{n_depots}};

}
parameters{ 
  
  positive_ordered[n_depots] TVKA;
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  positive_ordered[n_depots_with_delay] TVDUR_backward;
  
  vector[n_depots] TVFRAC;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  vector<lower = 0>[2] sigma;
  cholesky_factor_corr[2] L_Sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{
  
  vector[n_time_new] ipred;       // ipred at the new timepoints (no residual error)
  vector[n_time_new] pred;        // pred at the new timepoints
  vector[n_time_new] dv;          // dv at the new timepoints (with residual error)
  vector[n_subjects_new] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[n_subjects_new] c_max;   // Cmax between t1 and t2 (c_max_ss)
  vector[n_subjects_new] t_max;   // Tmax between t1 and t2, then subtract off t1
  vector[n_subjects_new] t_half;  // Thalf
  
  array[n_subjects_new] row_vector[n_depots] KA;
  vector[n_subjects_new] CL;
  vector[n_subjects_new] VC;
  array[n_subjects_new] row_vector[n_depots_with_delay] DUR;
  vector[n_subjects_new] KE;
  
  vector[n_depots_with_delay] TVDUR = reverse(TVDUR_backward);
 
  {
    row_vector[n_random] typical_values = 
                            append_col(append_col(TVKA', [TVCL, TVVC]), TVDUR');

    matrix[n_subjects_new, n_random] eta_new;
    matrix[n_subjects_new, n_random] theta_new;
    
    matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
    matrix[n_random, n_random] Omega = quad_form_diag(R, omega);
    
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
    
    real sigma_sq_p = Sigma[1, 1];
    real sigma_sq_a = Sigma[2, 2];
    real sigma_p_a = Sigma[1, 2];
    
    matrix[n_time_new, n_cmt] x_ipred;
    matrix[n_time_new, n_cmt] x_pred;
    
    array[n_time_new] real rate;
    array[n_time_new] real rate_p;
    
    array[n_cmt] real bioav = append_array(to_array_1d(TVFRAC), 
                                           rep_array(1.0, (n_cmt - n_depots))); 

    for(i in 1:n_subjects_new){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects_new) .* exp(eta_new));
    
    CL = col(theta_new, n_depots + 1);
    VC = col(theta_new, n_depots + 2);

    KE = CL ./ VC;

    for(j in 1:n_subjects_new){
      
      KA[j] = theta_new[j, 1:n_depots];
      DUR[j] = theta_new[j, (n_depots + 2 + 1):n_random];
      
      // The x_p values are here to make it easier to incorporate covariates later
      vector[n_depots] ka_p = TVKA;
      real cl_p = TVCL;
      real vc_p = TVVC;
      vector[n_depots_with_delay] dur_p = TVDUR;
                    
      for(i in subj_start[j]:subj_end[j]){

        if(cmt[i] <= n_depots_with_delay){
          
          rate[i] = amt[i]/DUR[j, cmt[i]];
          rate_p[i] = amt[i]/dur_p[cmt[i]];
          
        }else{
          rate[i] = 0;
          rate_p[i] = 0;
        }
        if(is_inf(rate[i])) rate[i] = 0;
        if(is_inf(rate_p[i])) rate_p[i] = 0;
      
      }
      
      x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(parallel_1cmt_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         append_array(to_array_1d(KA[j]), {CL[j], VC[j]}), 
                         bioav, tlag, x_r, x_i)';
                      
      ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], (n_depots + 1)] ./ VC[j];

      x_pred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(parallel_1cmt_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         append_array(to_array_1d(ka_p), {cl_p, vc_p}), 
                         bioav, tlag, x_r, x_i)';

      pred[subj_start[j]:subj_end[j]] = 
        x_pred[subj_start[j]:subj_end[j], (n_depots + 1)] ./ TVVC;
        
      auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], n_depots + 2]) / VC[j];
      c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], n_depots + 3]);
      t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], n_depots + 4]) - t_1;
      t_half[j] = log(2) ./ KE[j];
    }

    for(i in 1:n_time_new){
      if(ipred[i] == 0){
        dv[i] = 0;
      }else{
        real ipred_tmp = ipred[i];
        real sigma_tmp = sqrt(square(ipred_tmp) * sigma_sq_p + sigma_sq_a + 
                              2*ipred_tmp*sigma_p_a);
        dv[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
      }
    }
  }
}
