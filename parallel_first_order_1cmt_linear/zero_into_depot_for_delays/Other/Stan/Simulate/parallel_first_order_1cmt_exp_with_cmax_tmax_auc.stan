// One-compartment PK Model with Parallel absorption
// This is flexible so that any number of parallel absorption processes can take 
//   place (n_depots >= 2).
// A delay in absorption for each process is implemented with a zero-order 
//   distributive delay (like an infusion into the depot) to bring about a
//   delayed absorption. The fastest process may or may not have a delay (user's
//   choice). Control this with "delay_on_fastest_depot"
// There will be n_depots KA values for each subject. i.e., TVKA and KAi are 
//   vectors. KA is intended to be ordered so that 
//   KA_1 < KA_2 < ... < KA_n_depots. This doesn't do this exactly, but it will 
//   have TVKA_1 < TVKA_2 < ... < TVKA_n_depots and hope the individual effects 
//   don't mess that up
// There will be n_depots DUR values for each subject. i.e., TVDUR and DURi are 
//   vectors. DUR is intended to be ordered so that 
//   DUR_1 > DUR_2 > ... > DUR_n_depots. This doesn't do this exactly, but it 
//   will have TVDUR_1 > TVDUR_2 > ... > TVDUR_n_depots and hope the individual
//   effects don't mess that up
// FRAC is a simplex of length n_depots that tells how much of the dose goes into
//   each absorption process. There is no IIV on FRAC
// IIV on CL, VC, KA, DUR (full covariance matrix)
// General ODE solution using Torsten
// exponential error - DV = IPRED*exp(eps)
// Output includes individual Cmax over the whole time period, Tmax between t1 
//   and t2, AUC since 0 for every timepoint, and AUC between t1 and t2 (like a 
//   dosing interval)

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
  
  // An ODE that is written specifically for n_depots = 3. I don't ever use this
  //   one, but it is here to show what the flexible one above is doing.
  vector parallel_3_1cmt_ode(real t, vector y, array[] real params, 
                             array[] real x_r, array[] int x_i){
    
    real ka_1 = params[1];
    real ka_2 = params[2];
    real ka_3 = params[3];
    real cl = params[4];
    real vc = params[5];
    
    real ke = cl/vc;
    
    vector[4] dydt;

    dydt[1] = -ka_1*y[1];           // depot_slow
    dydt[2] = -ka_2*y[2];           // depot_medium
    dydt[3] = -ka_3*y[3];           // depot_fast
    dydt[4] = ka_1*y[1] + ka_2*y[2] + ka_3*y[3] - ke*y[4];  // central
    
    return dydt;
  }
  
}
data{
  
  int n_subjects;
  int n_total;                  
  
  array[n_total] real amt;
  array[n_total] int cmt;
  array[n_total] int evid;
  // array[n_total] real rate;
  array[n_total] real ii;
  array[n_total] int addl;
  array[n_total] int ss;
  array[n_total] real time;
  
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;
  
  int<lower = 2> n_depots; 
  int<lower = 0, upper = 1> delay_on_fastest_depot;
  
  positive_ordered[n_depots] TVKA;
  real<lower = 0> TVCL;
  real<lower = 0> TVVC;
  vector<lower = 0>[n_depots] TVDUR; // SHOULD BE DECREASING ORDERED
  simplex[n_depots] TVFRAC;
  
  vector<lower = 0>[n_depots] omega_ka;
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  vector<lower = 0>[n_depots] omega_dur;
  
  corr_matrix[2*n_depots + 2] R; // Correlation matrix before transforming to Omega.
                                 // Can in theory change this to having inputs for
                                 // cor_cl_vc, cor_cl_ka, ... and then construct the 
                                 // correlation matrix in transformed data, but it's easy
                                 // enough to do in R
  
  real<lower = 0> sigma;
  
  real<lower = 0> t_1;
  real<lower = t_1> t_2;
  
}
transformed data{
  
  int n_random = 2*n_depots + 2; // Number of random effects
  int n_cmt = n_depots + 4;      // Depot_1, ..., Depot_n, central, AUC_t1_t2, Cmax, Tmax
  
  vector<lower = 0>[n_depots] TVDUR_to_use = TVDUR;
  vector<lower = 0>[n_depots] omega_dur_to_use = omega_dur;
  
  if(delay_on_fastest_depot == 0){
    
    TVDUR_to_use[n_depots] = 0;
    omega_dur_to_use[n_depots] = 0;
    
  }
  
  vector[n_random] omega = append_row(
                             append_row(
                               append_row(to_vector(omega_ka), omega_cl), 
                               omega_vc), 
                             to_vector(omega_dur_to_use));
  
  matrix[n_random, n_random] L = cholesky_decompose(R);
  
  array[n_cmt] real bioav = append_array(to_array_1d(TVFRAC), 
                                           rep_array(1.0, (n_cmt - n_depots)));
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);  // Hardcoding. The delay comes from the infusion into the depots
  
  array[1, 2] real x_r = {{t_1, t_2}};
  
  array[1, 1] int x_i = {{n_depots}};
  
}
model{
  
}
generated quantities{
  
  vector[n_total] ipred; // concentration with no residual error
  vector[n_total] dv;    // concentration with residual error
  vector[n_subjects] auc_t1_t2;  // AUC from t1 up to t2 
  vector[n_subjects] c_max;      // Cmax between t1 and t2 
  vector[n_subjects] t_max;      // Tmax between t1 and t2, then subtract off t1
  vector[n_subjects] t_half;     // Thalf
  
  array[n_subjects] row_vector[n_depots] KA;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  array[n_subjects] row_vector[n_depots] DUR;
  vector[n_subjects] KE;
  
  {
  
    array[n_total] real rate;
    
    row_vector[n_random] typical_values = 
                            append_col(append_col(TVKA', [TVCL, TVVC]), TVDUR_to_use');
    
    matrix[n_subjects, n_random] eta;   
    matrix[n_subjects, n_random] theta; 
  
    matrix[n_total, n_cmt] x_ipred;
    
    for(i in 1:n_subjects){
      eta[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                           diag_pre_multiply(omega, L))';
    }
    theta = (rep_matrix(typical_values, n_subjects) .* exp(eta));

    CL = col(theta, n_depots + 1);
    VC = col(theta, n_depots + 2);
    KE = CL ./ VC;
    
    for(j in 1:n_subjects){
      
      KA[j] = theta[j, 1:n_depots];
      DUR[j] = theta[j, (n_depots + 2 + 1):n_random];
      
      for(i in subj_start[j]:subj_end[j]){
        
        if(cmt[i] <= n_depots){
          rate[i] = amt[i]/DUR[j, cmt[i]];
        }else{
          rate[i] = 0;
        }
        if(is_inf(rate[i])) rate[i] = 0;
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
                    
      auc_t1_t2[j] = max(x_ipred[subj_start[j]:subj_end[j], n_depots + 2]) / VC[j];
      c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], n_depots + 3]);
      t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], n_depots + 4]) - t_1;
      t_half[j] = log(2) ./ KE[j];
    
    }

    for(i in 1:n_total){
      if(ipred[i] == 0){
         dv[i] = 0;
      }else{
        dv[i] = lognormal_rng(log(ipred[i]), sigma);
      }
    }
  }
}


