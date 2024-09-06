// IV infusion
// Three-compartment PK Model
// IIV on CL, VC, Q1, VP1, Q2, and VP2 (full covariance matrix)
// proportional plus additive error - DV = IPRED*(1 + eps_p) + eps_a
// General ODE solution using Torsten to get out individual estimates of AUC. To
//   get Cmax and Tmax, make sure you simulate at the end of the infusion
// Predictions are generated from a normal that is truncated below at 0

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  vector iv_3cmt_ode(real t, vector y, array[] real params, 
                     array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real q1 = params[3];
    real vp1 = params[4];
    real q2 = params[5];
    real vp2 = params[6];
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real ke = cl/vc;
    real k_cp1 = q1/vc;
    real k_p1c = q1/vp1;
    real k_cp2 = q2/vc;
    real k_p2c = q2/vp2;
    
    vector[5] dydt;

    dydt[1] = -(ke + k_cp1 + k_cp2)*y[1] + k_p1c*y[2] + k_p2c*y[3];  // central
    dydt[2] = k_cp1*y[1] - k_p1c*y[2];                               // peripheral 1
    dydt[3] = k_cp2*y[1] - k_p2c*y[3];                               // peripheral 2
    dydt[4] = y[1];                                                  // AUC
    dydt[5] = t >= t_1 && t <= t_2 ? y[1] : 0;                       // AUC_t_1-t_2
    
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
  array[n_time_new] real rate;
  array[n_time_new] real ii;
  array[n_time_new] int addl;
  array[n_time_new] int ss;
  array[n_subjects_new] int subj_start;
  array[n_subjects_new] int subj_end;
  
  real<lower = 0> t_1;   // Time at which to start AUC calculations (AUC_t1_t2, ...)
  real<lower = t_1> t_2; // Time at which to end AUC calculations (AUC_t1_t2, ...)
 
}
transformed data{ 
  
  int n_random = 6;    // Number of random effects
  int n_cmt = 3;       // Number of ODEs in PK model (central, peripheral1, peripheral2)
  int n_cmt_extra = 2; // Number of ODEs for AUCs
  
  array[n_cmt + n_cmt_extra] real bioav = rep_array(1.0, n_cmt + n_cmt_extra);
  array[n_cmt + n_cmt_extra] real tlag = rep_array(0.0, n_cmt + n_cmt_extra);
  
  array[1, 2] real x_r = {{t_1, t_2}};

}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVQ1;       
  real<lower = 0> TVVP1;
  real<lower = 0> TVQ2;       
  real<lower = 0> TVVP2;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  vector<lower = 0>[2] sigma;
  cholesky_factor_corr[2] L_Sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{
  
  vector[n_time_new] ipred;      // ipred for the observed individuals at the new timepoints
  vector[n_time_new] pred;       // pred for the observed individuals at the new timepoints
  vector[n_time_new] dv;         // dv for the observed individuals at the new timepoints
  vector[n_time_new] auc;        // AUC for the observed individuals at the new timepoints
  vector[n_subjects_new] auc_t1_t2;  // AUC from t1 up to t2 (AUC_ss)
  vector[n_subjects_new] t_half_alpha;    // alpha half-life
  vector[n_subjects_new] t_half_beta;     // beta half-life
  vector[n_subjects_new] t_half_terminal; // terminal half-life
  
  vector[n_subjects_new] CL;
  vector[n_subjects_new] VC;
  vector[n_subjects_new] Q1;
  vector[n_subjects_new] VP1;
  vector[n_subjects_new] Q2;
  vector[n_subjects_new] VP2;
 
  {
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ1, TVVP1,
                                                                     TVQ2, TVVP2});
                                                                     
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
    
    real sigma_sq_p = Sigma[1, 1];
    real sigma_sq_a = Sigma[2, 2];
    real sigma_p_a = Sigma[1, 2];
    
    matrix[n_subjects_new, n_random] eta_new;
    matrix[n_subjects_new, n_random] theta_new;

    matrix[n_time_new, n_cmt + n_cmt_extra] x_pred;
    matrix[n_time_new, n_cmt + n_cmt_extra] x_ipred;
    
    for(i in 1:n_subjects_new){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects_new) .* exp(eta_new));
    
    CL = col(theta_new, 1);
    VC = col(theta_new, 2);
    Q1 = col(theta_new, 3);
    VP1 = col(theta_new, 4);
    Q2 = col(theta_new, 5);
    VP2 = col(theta_new, 6);

    for(j in 1:n_subjects_new){
      
      real cl_p = TVCL;   // * wt_adjustment_cl * egfr_adjustment_cl;
      real vc_p = TVVC;   // * wt_adjustment_vc * race_asian_adjustment_vc;
      real q1_p = TVQ1;   // * wt_adjustment_q1;
      real vp1_p = TVVP1; // * wt_adjustment_vp1;
      real q2_p = TVQ2;   // * wt_adjustment_q2;
      real vp2_p = TVVP2; // * wt_adjustment_vp2;
      
      real ke = CL[j]/VC[j];
      real k_cp1 = Q1[j]/VC[j];
      real k_p1c = Q1[j]/VP1[j];
      real k_cp2 = Q2[j]/VC[j];
      real k_p2c = Q2[j]/VP2[j];
      
      real jay = k_cp1 + ke + k_p1c + k_p2c + k_cp2;
      real kay = k_cp1*k_p2c + ke*k_p1c + ke*k_p2c + k_p1c*k_p2c + k_cp2*k_p1c;
      real ell = ke*k_p1c*k_p2c;
      real m = (3*kay - square(jay))/3;
      real n = (2*jay^3 - 9*jay*kay + 27*ell)/27;
      real Q = square(n)/4 + m^3/27;
      real alpha = sqrt(-Q);
      real beta = -n/2;
      real rho = hypot(beta, alpha);  // sqrt(square(beta) + square(alpha))
      real delta = atan2(alpha, beta);
      
      real lambda_1 = jay/3 + cbrt(rho)*(cos(delta/3) + sqrt(3)*sin(delta/3));
      real lambda_2 = jay/3 + cbrt(rho)*(cos(delta/3) - sqrt(3)*sin(delta/3));
      real lambda_3 = jay/3 - 2*cbrt(rho)*cos(delta/3);
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(iv_3cmt_ode,
                       n_cmt + n_cmt_extra,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {CL[j], VC[j], Q1[j], VP1[j], Q2[j], VP2[j]}, 
                       bioav, tlag, x_r)';
                      
      ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], 1] ./ VC[j];

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(iv_3cmt_ode,
                       n_cmt + n_cmt_extra,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {cl_p, vc_p, q1_p, vp1_p, q2_p, vp2_p}, 
                       bioav, tlag, x_r)';
                      

      pred[subj_start[j]:subj_end[j]] = 
        x_pred[subj_start[j]:subj_end[j], 1] ./ vc_p;
      
      auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 4] ./ VC[j];
      
      auc_t1_t2[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]) / VC[j];

      t_half_alpha[j] = log(2)/lambda_1;
      t_half_beta[j] = log(2)/lambda_2;
      t_half_terminal[j] = log(2)/lambda_3;
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

