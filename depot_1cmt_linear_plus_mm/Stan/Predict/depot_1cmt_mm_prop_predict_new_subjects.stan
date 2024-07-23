// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// IIV on VC, VMAX, KM, KA
// proportional error - DV = IPRED*(1 + eps_p)
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
  
  vector depot_1cmt_mm_ode(real t, vector y, array[] real params, 
                           array[] real x_r, array[] int x_i){
  
    real vc = params[1];
    real vmax = params[2];
    real km = params[3];
    real ka = params[4];
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real conc = y[2]/vc;
    
    real slope = ka*y[1] - vmax*conc/(km + conc);
    real x = slope > 0 && conc > y[5] ? slope/vc : 0;
    real z = t <= t_1 || (slope > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    vector[6] dydt;
    
    dydt[1] = -ka*y[1];
    dydt[2] = slope;
    dydt[3] = y[2];                                // AUC
    dydt[4] = t >= t_1 && t <= t_2 ? y[2] : 0;     // AUC_t_1-t_2
    dydt[5] = x;                                   // C_max
    dydt[6] = z;                                   // t_max
    
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
  
  real<lower = 0> t_1;   // Time at which to start SS calculations (AUC_ss, C_max_ss, ...)
  real<lower = t_1> t_2; // Time at which to end SS calculations (AUC_ss, C_max_ss, ...)
 
}
transformed data{ 
  
  int n_random = 4; // Number of random effects
  int n_cmt = 2;
  int n_cmt_extra = 4; // AUC, Cmax, Tmax compartments
  
  array[n_cmt + n_cmt_extra] real bioav = rep_array(1.0, n_cmt + n_cmt_extra);
  array[n_cmt + n_cmt_extra] real tlag = rep_array(0.0, n_cmt + n_cmt_extra);
  
  array[1, 2] real x_r = {{t_1, t_2}};

}
parameters{ 
  
  real<lower = 0> TVVC;
  real<lower = 0> TVVMAX;       
  real<lower = 0> TVKM;
  real<lower = 0> TVKA;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{
  
  vector[n_time_new] ipred;       // ipred for the new individuals at the new timepoints
  vector[n_time_new] pred;        // pred for the new individuals at the new timepoints
  vector[n_time_new] dv;          // dv for the new individuals at the new timepoints
  vector[n_time_new] auc;         // auc for the new individuals from time 0 to the new timepoint
  vector[n_subjects_new] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[n_subjects_new] c_max;   // Cmax between t1 and t2 (c_max_ss)
  vector[n_subjects_new] t_max;   // Tmax between t1 and t2, then subtract off t1
  
  vector[n_subjects_new] VC;
  vector[n_subjects_new] VMAX;
  vector[n_subjects_new] KM;
  vector[n_subjects_new] KA;

 
  {
    
    row_vector[n_random] typical_values = to_row_vector({TVVC, TVVMAX, TVKM, TVKA});
    
    matrix[n_subjects_new, n_random] eta_new;
    matrix[n_subjects_new, n_random] theta_new;
    
    matrix[n_time_new, n_cmt + n_cmt_extra] x_pred;
    matrix[n_time_new, n_cmt + n_cmt_extra] x_ipred;
    
    for(i in 1:n_subjects_new){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects_new) .* exp(eta_new));
    
    VC = col(theta_new, 1);
    VMAX = col(theta_new, 2);
    KM = col(theta_new, 3);
    KA = col(theta_new, 4);

    for(j in 1:n_subjects_new){
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_1cmt_mm_ode,
                       n_cmt + n_cmt_extra,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {VC[j], VMAX[j], KM[j], KA[j]},
                       bioav, tlag, x_r)';

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_1cmt_mm_ode,
                       n_cmt + n_cmt_extra,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {TVVC, TVVMAX, TVKM, TVKA},
                       bioav, tlag, x_r)';
                      
      ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];

      pred[subj_start[j]:subj_end[j]] = 
        x_pred[subj_start[j]:subj_end[j], 2] ./ TVVC;
        
      auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 3] ./ VC[j];
        
      auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], 4]) / VC[j];
      c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]);
      t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 6]) - t_1;

    }

  
    for(i in 1:n_time_new){
      if(ipred[i] == 0){
        dv[i] = 0;
      }else{
        real ipred_tmp = ipred[i];
        real sigma_tmp = ipred_tmp*sigma_p;
        dv[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
      }
    }
  }
}
