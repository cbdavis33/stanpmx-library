// IV infusion
// Two-compartment PK Model
// IIV on CL, VC, Q, and VP (full covariance matrix)
// proportional plus additive error - DV = IPRED*(1 + eps_p) + eps_a
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
  int n_time_new;
  array[n_time_new] real time;
  array[n_time_new] real amt;
  array[n_time_new] int cmt;
  array[n_time_new] int evid;
  array[n_time_new] real rate;
  array[n_time_new] real ii;
  array[n_time_new] int addl;
  array[n_time_new] int ss;
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;
  
  real<lower = 0> t_1;   // Time at which to start SS calculations (AUC_ss, ...)
  real<lower = t_1> t_2; // Time at which to end SS calculations (AUC_ss, ...)
  
  int<lower = 0, upper = 1> want_auc_cmax; // Want AUC and Cmax? If so, it'll 
                                           // use the ODE solution. Otherwise,
                                           // it'll use the matrix-exponential 
                                           // solution (and be faster)
 
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
  
  vector<lower = 0>[2] sigma;
  cholesky_factor_corr[2] L_Sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{
  
  vector[n_time_new] pred;       // f(TVs, x, eta = 0) 
  vector[n_time_new] epred_stan; // f(TVs, x, eta = eta_new), eta_new ~ multi_normal(0, Omega) 
  vector[n_time_new] epred;      // epred_stan + error
  vector[n_time_new] ipred;      // f(TVs, x, eta = eta_i), eta_i are etas for observed subjects
  vector[n_time_new] dv;         // ipred + error
  
  vector[want_auc_cmax ? n_subjects : 0] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[want_auc_cmax ? n_time_new : 0] auc;     // AUC for the observed individuals at the new timepoints
  vector[n_subjects] t_half_alpha;                // alpha half-life
  vector[n_subjects] t_half_terminal;             // terminal half-life
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;
 
  {
    
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP});
    
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
    
    real sigma_sq_p = Sigma[1, 1];
    real sigma_sq_a = Sigma[2, 2];
    real sigma_p_a = Sigma[1, 2]; 

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));

    matrix[n_subjects, n_random] eta_new;
    matrix[n_subjects, n_random] theta_new;
    
    for(i in 1:n_subjects){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects) .* exp(eta_new));
    
    vector[n_subjects] CL_new;
    vector[n_subjects] VC_new;
    vector[n_subjects] Q_new;
    vector[n_subjects] VP_new;
    
    matrix[n_time_new, 2] x_pred;
    matrix[n_time_new, 2] x_epred;
    matrix[n_time_new, n_cmt] x_ipred;
    
    vector[n_subjects] alpha;
    vector[n_subjects] beta;
    
    for(j in 1:n_subjects){
      
      row_vector[n_random] theta_j = theta[j]; // access the parameters for subject j
      row_vector[n_random] theta_j_new = theta_new[j]; // access the parameters for subject j's epred
      
      real cl_p = TVCL;
      real vc_p = TVVC;
      real q_p = TVQ;
      real vp_p = TVVP;
      
      CL[j] = theta_j[1];
      VC[j] = theta_j[2];
      Q[j] = theta_j[3];
      VP[j] = theta_j[4];
      
      CL_new[j] = theta_j_new[1];
      VC_new[j] = theta_j_new[2];
      Q_new[j] = theta_j_new[3];
      VP_new[j] = theta_j_new[4];
      
      matrix[2, 2] K_epred = rep_matrix(0, 2, 2);
      matrix[2, 2] K_p = rep_matrix(0, 2, 2);
      
      real ke_new = CL_new[j]/VC_new[j];
      real k_cp_new = Q_new[j]/VC_new[j];
      real k_pc_new = Q_new[j]/VP_new[j];
    
      K_epred[1, 1] = -(ke_new + k_cp_new);
      K_epred[1, 2] = k_pc_new;
      K_epred[2, 1] = k_cp_new;
      K_epred[2, 2] = -k_pc_new;
      
      x_epred[subj_start[j]:subj_end[j], ] =
        pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         K_epred, bioav, tlag)';
                         
      epred_stan[subj_start[j]:subj_end[j]] = 
        x_epred[subj_start[j]:subj_end[j], 1] ./ VC_new[j];
      
      real ke_p = cl_p/vc_p;
      real k_cp_p = q_p/vc_p;
      real k_pc_p = q_p/vp_p;
    
      K_p[1, 1] = -(ke_p + k_cp_p);
      K_p[1, 2] = k_pc_p;
      K_p[2, 1] = k_cp_p;
      K_p[2, 2] = -k_pc_p;
      
      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         K_p, bioav, tlag)';

      pred[subj_start[j]:subj_end[j]] = 
        x_pred[subj_start[j]:subj_end[j], 1] ./ vc_p;
      
      if(want_auc_cmax == 1){
        
        x_ipred[subj_start[j]:subj_end[j],] =
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
                        
        ipred[subj_start[j]:subj_end[j]] = 
          x_ipred[subj_start[j]:subj_end[j], 1] ./ VC[j];
          
        auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], 4]) / VC[j];
        auc[subj_start[j]:subj_end[j]] =
                                x_ipred[subj_start[j]:subj_end[j], 3] ./ VC[j];
          
      }else{
        
        matrix[2, 2] K = rep_matrix(0, 2, 2);
        
        real ke = CL[j]/VC[j];
        real k_cp = Q[j]/VC[j];
        real k_pc = Q[j]/VP[j];
    
        K[1, 1] = -(ke + k_cp);
        K[1, 2] = k_pc;
        K[2, 1] = k_cp;
        K[2, 2] = -k_pc;
      
        x_ipred[subj_start[j]:subj_end[j], ] =
          pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           K, bioav, tlag)';

        ipred[subj_start[j]:subj_end[j]] = 
          x_ipred[subj_start[j]:subj_end[j], 1] ./ VC[j];
        
      }

    }
    
    alpha = 0.5*(CL./VC + Q./VC + Q./VP + 
                   sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
    beta = 0.5*(CL./VC + Q./VC + Q./VP - 
                   sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
                
    t_half_alpha = log(2)/alpha;
    t_half_terminal = log(2)/beta;
  
    for(i in 1:n_time_new){
      
      if(epred_stan[i] == 0){
        epred[i] = 0;
      }else{
        real epred_tmp = epred_stan[i];
        real sigma_tmp_e = sqrt(square(epred_tmp) * sigma_sq_p + sigma_sq_a + 
                                2*epred_tmp*sigma_p_a);
        epred[i] = normal_lb_rng(epred_tmp, sigma_tmp_e, 0.0);
      }
      
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

