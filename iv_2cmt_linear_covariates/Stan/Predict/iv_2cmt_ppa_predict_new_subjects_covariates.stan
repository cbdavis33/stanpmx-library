// IV infusion
// Two-compartment PK Model
// IIV on CL, VC, Q, and VP (full covariance matrix)
// proportional plus additive error - DV = IPRED*(1 + eps_p) + eps_a
// General ODE solution using Torsten to get out individual estimates of AUC. To
//   get Cmax and Tmax, make sure you simulate at the end of the infusion
// Predictions are generated from a normal that is truncated below at 0
// Covariates: 
//   1) Body Weight on CL, VC, Q, VP - (wt/70)^theta
//   2) Race = Asian on VC (0/1) - exp(theta*race_asian)
//   3) eGFR on CL (continuous) - (eGFR/90)^theta

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
  
  array[n_subjects_new] real<lower = 0> wt;                   // baseline bodyweight (kg)
  array[n_subjects_new] int<lower = 0, upper = 1> race_asian; // race = Asian
  array[n_subjects_new] real<lower = 0> egfr;                 // eGFR
  
  real<lower = 0> t_1;   // Time at which to start SS calculations (AUC_ss, ...)
  real<lower = t_1> t_2; // Time at which to end SS calculations (AUC_ss, ...)
  
}
transformed data{ 
  
  int n_random = 4; // Number of random effects
  int n_cmt = 4;    // Number of compartments (depot, central, AUC, AUC_ss)
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt);
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  array[1, 2] real x_r = {{t_1, t_2}};

}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVQ;       
  real<lower = 0> TVVP;
  
  real theta_cl_wt;         // allometric scaling coefficient for wt on clearance
  real theta_vc_wt;         // allometric scaling coefficient for wt on VC
  real theta_q_wt;          // allometric scaling coefficient for wt on Q
  real theta_vp_wt;         // allometric scaling coefficient for wt on VP
  real theta_vc_race_asian; // effect of race (asian) on VC 
  real theta_cl_egfr;       // effect of eGFR on clearance
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  vector<lower = 0>[2] sigma;
  cholesky_factor_corr[2] L_Sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{

  vector[n_time_new] ipred;   // ipred for the observed individuals at the new timepoints
  vector[n_time_new] pred;    // pred for the observed individuals at the new timepoints
  vector[n_time_new] dv;      // dv for the observed individuals at the new timepoints
  vector[n_time_new] auc;     // AUC for the observed individuals at the new timepoints
  vector[n_subjects_new] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[n_subjects_new] t_half_alpha;    // alpha half-life
  vector[n_subjects_new] t_half_terminal; // terminal half-life
  
  vector[n_subjects_new] CL;
  vector[n_subjects_new] VC;
  vector[n_subjects_new] Q;
  vector[n_subjects_new] VP;

  {
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP});
    
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
    
    real sigma_sq_p = Sigma[1, 1];
    real sigma_sq_a = Sigma[2, 2];
    real sigma_p_a = Sigma[1, 2];
    
    matrix[n_subjects_new, n_random] eta_new;
    matrix[n_subjects_new, n_random] theta_new;
    matrix[n_time_new, n_cmt] x_ipred;
    matrix[n_time_new, n_cmt] x_pred;
    
    vector[n_subjects_new] alpha;
    vector[n_subjects_new] beta;

    for(i in 1:n_subjects_new){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects_new) .* exp(eta_new));

    for(j in 1:n_subjects_new){
      
      real wt_over_70 = wt[j]/70;
      real wt_adjustment_cl = wt_over_70^theta_cl_wt;
      real wt_adjustment_vc = wt_over_70^theta_vc_wt;
      real wt_adjustment_q = wt_over_70^theta_q_wt;
      real wt_adjustment_vp = wt_over_70^theta_vp_wt;
      real race_asian_adjustment_vc = exp(theta_vc_race_asian*race_asian[j]);
      real egfr_adjustment_cl = (egfr[j]/90)^theta_cl_egfr;
      
      real cl_p = TVCL * wt_adjustment_cl * egfr_adjustment_cl;
      real vc_p = TVVC * wt_adjustment_vc * race_asian_adjustment_vc;
      real q_p = TVQ * wt_adjustment_q;
      real vp_p = TVVP * wt_adjustment_vp;
      
      row_vector[n_random] theta_j = theta_new[j]; // access the parameters for subject j
      CL[j] = theta_j[1] * wt_adjustment_cl * egfr_adjustment_cl;
      VC[j] = theta_j[2] * wt_adjustment_vc * race_asian_adjustment_vc;
      Q[j] = theta_j[3] * wt_adjustment_q;
      VP[j] = theta_j[4] * wt_adjustment_vp;
      
      alpha[j] = 0.5*(CL[j]/VC[j] + Q[j]/VC[j] + Q[j]/VP[j] + 
                      sqrt((CL[j]/VC[j] + Q[j]/VC[j] + Q[j]/VP[j])^2 - 
                      4*CL[j]/VC[j]*Q[j]/VP[j]));
      beta[j] = 0.5*(CL[j]/VC[j] + Q[j]/VC[j] + Q[j]/VP[j] - 
                     sqrt((CL[j]/VC[j] + Q[j]/VC[j] + Q[j]/VP[j])^2 - 
                     4*CL[j]/VC[j]*Q[j]/VP[j]));
      
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
                       {CL[j], VC[j], Q[j], VP[j]}, bioav, tlag, x_r)';
                      
      ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], 1] ./ VC[j];

      x_pred[subj_start[j]:subj_end[j],] =
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
                       {cl_p, vc_p, q_p, vp_p}, bioav, tlag, x_r)';

      pred[subj_start[j]:subj_end[j]] = 
        x_pred[subj_start[j]:subj_end[j], 1] ./ vc_p;
        
      auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 3] ./ VC[j];
      
      auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], 4]) / VC[j];

      t_half_alpha[j] = log(2)/alpha[j];
      t_half_terminal[j] = log(2)/beta[j];
    
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


