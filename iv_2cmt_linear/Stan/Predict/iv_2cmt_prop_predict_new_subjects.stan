// IV infusion
// Two-compartment PK Model
// IIV on CL, VC, Q, and VP (full covariance matrix)
// proportional error - DV = IPRED*(1 + eps_p)
// General ODE solution using Torsten will get out individual estimates of AUC.
//   Matrix-exponential will not. To get Cmax and Tmax, make sure you simulate 
//   at the end of the infusion
// Predictions are generated from a normal that is truncated below at 0

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  vector iv_2cmt_ode(real t, vector y, array[] real params, 
                     array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real q = params[3];
    real vp = params[4];
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;
    
    vector[4] dydt;

    dydt[1] = -(ke + k_cp)*y[1] + k_pc*y[2];   // central
    dydt[2] = k_cp*y[1] - k_pc*y[2];           // peripheral
    dydt[3] = y[1];                            // AUC
    dydt[4] = t >= t_1 && t <= t_2 ? y[1] : 0; // AUC_t_1-t_2
    
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
  
  int<lower = 0, upper = 1> want_auc_cmax; // Want AUC and Cmax? If so, it'll 
                                           // use the ODE solution. Otherwise,
                                           // it'll use the analytical solution (and be faster)
  
}
transformed data{ 
  
  int n_random = 4; // Number of random effects
  int n_cmt = want_auc_cmax ? 4 : 2; // Number of compartments - central, peripheral (AUC, AUC_ss)
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt);
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  array[1, 2] real x_r = {{t_1, t_2}};

}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVQ;       
  real<lower = 0> TVVP;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{

  vector[n_time_new] epred_stan; // f(TVs, x, eta = eta_new), eta_new ~ multi_normal(0, Omega) 
  vector[n_time_new] epred;      // epred_stan + error
  
  vector[want_auc_cmax ? n_subjects_new : 0] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[want_auc_cmax ? n_time_new : 0] auc;         // AUC for the new individuals at the new timepoints
  vector[n_subjects_new] t_half_alpha;                // alpha half-life
  vector[n_subjects_new] t_half_terminal;             // terminal half-life
  
  vector[n_subjects_new] CL;
  vector[n_subjects_new] VC;
  vector[n_subjects_new] Q;
  vector[n_subjects_new] VP;

  {
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP});
    
    matrix[n_subjects_new, n_random] eta_new;
    matrix[n_subjects_new, n_random] theta_new;
    matrix[n_time_new, n_cmt] x_epred;

    for(i in 1:n_subjects_new){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects_new) .* exp(eta_new));
    
    vector[n_subjects_new] alpha;
    vector[n_subjects_new] beta;

    for(j in 1:n_subjects_new){
      
      row_vector[n_random] theta_j_new = theta_new[j]; // access the parameters for subject j's epred
      
      CL[j] = theta_j_new[1];
      VC[j] = theta_j_new[2];
      Q[j] = theta_j_new[3];
      VP[j] = theta_j_new[4];
    
      if(want_auc_cmax == 1){
        
        x_epred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(iv_2cmt_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j], Q[j], VP[j]}, 
                         bioav, tlag, x_r)';
                        
        epred_stan[subj_start[j]:subj_end[j]] = 
          x_epred[subj_start[j]:subj_end[j], 1] ./ VC[j];
          
        auc_ss[j] = max(x_epred[subj_start[j]:subj_end[j], 4]) / VC[j];
        auc[subj_start[j]:subj_end[j]] =
                                x_epred[subj_start[j]:subj_end[j], 3] ./ VC[j];
          
      }else{
        
        matrix[2, 2] K = rep_matrix(0, 2, 2);
        
        real ke = CL[j]/VC[j];
        real k_cp = Q[j]/VC[j];
        real k_pc = Q[j]/VP[j];
    
        K[1, 1] = -(ke + k_cp);
        K[1, 2] = k_pc;
        K[2, 1] = k_cp;
        K[2, 2] = -k_pc;
      
        x_epred[subj_start[j]:subj_end[j], ] =
          pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           K, bioav, tlag)';
        
        epred_stan[subj_start[j]:subj_end[j]] = 
          x_epred[subj_start[j]:subj_end[j], 1] ./ VC[j];
        
      }
  
      alpha = 0.5*(CL./VC + Q./VC + Q./VP + 
                   sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
      beta = 0.5*(CL./VC + Q./VC + Q./VP - 
                  sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
                  
      t_half_alpha = log(2)/alpha;
      t_half_terminal = log(2)/beta;
    
    }

    for(i in 1:n_time_new){
      if(epred_stan[i] == 0){
        epred[i] = 0;
      }else{
        real epred_tmp = epred_stan[i];
        real sigma_tmp_e = epred_tmp*sigma_p;
        epred[i] = normal_lb_rng(epred_tmp, sigma_tmp_e, 0.0);
      }
    }
  }
}


