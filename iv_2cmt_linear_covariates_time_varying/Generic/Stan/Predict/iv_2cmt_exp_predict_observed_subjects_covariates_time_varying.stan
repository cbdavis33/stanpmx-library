// IV infusion
// Two-compartment PK Model
// IIV on CL, VC, Q, and VP (full covariance matrix)
// exponential error - DV = IPRED*exp(eps)
// General ODE solution using Torsten to get out individual estimates of AUC, 
//   Cmax, Tmax, ... is an option, but the user can choose to only calculate 
//   IPRED and DV
// Covariates - this file is generic, so all can be time-varying or constant. 
//   The key will be that the input for each covariate is length n_ttime_new and 
//   not of length n_subjects (length = n_subjects implies that covariate is 
//   constant): 
//   1) Body Weight on CL, VC, Q, VP - (wt/70)^theta
//   2) Race = Asian on VC (0/1) - exp(theta*race_asian) (obviously race is 
//        constant, but this is for an example. It'll just be constant for 
//        length n_total)
//   3) eGFR on CL (continuous) - (eGFR/90)^theta

functions{
  
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

    dydt[1] = -(ke + k_cp)*y[1] + k_pc*y[2];      // central
    dydt[2] = k_cp*y[1] - k_pc*y[2];              // peripheral
    dydt[3] = y[1]/vc;                            // AUC
    dydt[4] = t >= t_1 && t <= t_2 ? y[1]/vc : 0; // AUC_t_1-t_2
    
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
  
  vector<lower = 0>[n_time_new] wt;                     // bodyweight (kg)
  vector<lower = 0, upper = 1>[n_time_new] race_asian;  // race = Asian
  vector<lower = 0>[n_time_new] egfr;                   // eGFR
  
  real<lower = 0> t_1;   // Time at which to start SS calculations (AUC_ss, ...)
  real<lower = t_1> t_2; // Time at which to end SS calculations (AUC_ss, ...)
  
  int <lower = 0, upper = 0> want_auc_cmax; // 0 => only calculate concentrations
                                            // 1 => concentrations and auc, ...
                                            // For now, the auc 
                                            // calculations aren't working as expected
                                            // with the time-varying parameters,
                                            // so this must be 0
 
}
transformed data{ 
  
  int n_random = 4; // Number of random effects
  int n_cmt = want_auc_cmax ? 4 : 3; // Number of compartments - central, peripheral (AUC, AUC_ss))
  
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
  
  real<lower = 0> sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{
  
  vector[n_time_new] ipred;   // ipred for the observed individuals at the new timepoints
  vector[n_time_new] pred;    // pred for the observed individuals at the new timepoints
  vector[n_time_new] dv;      // dv for the observed individuals at the new timepoints
  vector[want_auc_cmax ? n_subjects : 0] auc_ss;  // AUC from t1 up to t2 (AUC_ss)

  vector[n_time_new] CL;
  vector[n_time_new] VC;
  vector[n_time_new] Q;
  vector[n_time_new] VP;
 
  {
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP});

    matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
    matrix[n_random, n_random] Omega = quad_form_diag(R, omega);

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));

    matrix[n_time_new, 3] x_pred;
    matrix[n_time_new, n_cmt] x_ipred;
    
    vector[n_time_new] cl_p;
    vector[n_time_new] vc_p;
    vector[n_time_new] q_p;
    vector[n_time_new] vp_p;

    for(j in 1:n_subjects){
      
      int n_total_subj = subj_end[j] - subj_start[j] + 1;
      vector[n_total_subj] wt_over_70 = wt[subj_start[j]:subj_end[j]] ./ 70;
      vector[n_total_subj] wt_adjustment_cl = pow(wt_over_70, theta_cl_wt);
      vector[n_total_subj] wt_adjustment_vc = pow(wt_over_70, theta_vc_wt);
      vector[n_total_subj] wt_adjustment_q = pow(wt_over_70, theta_q_wt);
      vector[n_total_subj] wt_adjustment_vp = pow(wt_over_70, theta_vp_wt);
      vector[n_total_subj] race_asian_adjustment_vc = 
                 exp(theta_vc_race_asian*race_asian[subj_start[j]:subj_end[j]]);
      vector[n_total_subj] egfr_adjustment_cl = 
                      pow(egfr[subj_start[j]:subj_end[j]] ./ 90, theta_cl_egfr);
                      
      array[n_total_subj, want_auc_cmax ? n_random : n_random + 1] real params_to_input;
      array[n_total_subj, n_random + 1] real params_to_input_p;
                      
      row_vector[n_random] theta_j = theta[j]; // access the parameters for subject j
      CL[subj_start[j]:subj_end[j]] = theta_j[1] .* wt_adjustment_cl .* egfr_adjustment_cl;
      VC[subj_start[j]:subj_end[j]] = theta_j[2] .* wt_adjustment_vc .* race_asian_adjustment_vc;
      Q[subj_start[j]:subj_end[j]] = theta_j[3] * wt_adjustment_q;
      VP[subj_start[j]:subj_end[j]] = theta_j[4] * wt_adjustment_vp;
      
      cl_p[subj_start[j]:subj_end[j]] = TVCL .* wt_adjustment_cl .* egfr_adjustment_cl;
      vc_p[subj_start[j]:subj_end[j]] = TVVC .* wt_adjustment_vc .* race_asian_adjustment_vc;
      q_p[subj_start[j]:subj_end[j]] = TVQ * wt_adjustment_q;
      vp_p[subj_start[j]:subj_end[j]] = TVVP * wt_adjustment_vp;
      
      params_to_input_p[, 1] = to_array_1d(cl_p[subj_start[j]:subj_end[j]]);
      params_to_input_p[, 2] = to_array_1d(q_p[subj_start[j]:subj_end[j]]);
      params_to_input_p[, 3] = to_array_1d(vc_p[subj_start[j]:subj_end[j]]);
      params_to_input_p[, 4] = to_array_1d(vp_p[subj_start[j]:subj_end[j]]);
      params_to_input_p[, 5] = zeros_array(n_total_subj);
      
      if(want_auc_cmax == 1){
        
        params_to_input[, 1] = to_array_1d(CL[subj_start[j]:subj_end[j]]);
        params_to_input[, 2] = to_array_1d(VC[subj_start[j]:subj_end[j]]);
        params_to_input[, 3] = to_array_1d(Q[subj_start[j]:subj_end[j]]);
        params_to_input[, 4] = to_array_1d(VP[subj_start[j]:subj_end[j]]);
      
        params_to_input_p[, 1] = to_array_1d(cl_p[subj_start[j]:subj_end[j]]);
        params_to_input_p[, 2] = to_array_1d(vc_p[subj_start[j]:subj_end[j]]);
        params_to_input_p[, 3] = to_array_1d(q_p[subj_start[j]:subj_end[j]]);
        params_to_input_p[, 4] = to_array_1d(vp_p[subj_start[j]:subj_end[j]]);
        
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
                         params_to_input, bioav, tlag, x_r)';
                         
        ipred[subj_start[j]:subj_end[j]] = 
          x_ipred[subj_start[j]:subj_end[j], 1] ./ VC[subj_start[j]:subj_end[j]];
                         
        auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], 4]);
                         
      }else{
        
        params_to_input[, 1] = to_array_1d(CL[subj_start[j]:subj_end[j]]);
        params_to_input[, 2] = to_array_1d(Q[subj_start[j]:subj_end[j]]);
        params_to_input[, 3] = to_array_1d(VC[subj_start[j]:subj_end[j]]);
        params_to_input[, 4] = to_array_1d(VP[subj_start[j]:subj_end[j]]);
        params_to_input[, 5] = zeros_array(n_total_subj);
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           params_to_input)';
                           
        ipred[subj_start[j]:subj_end[j]] = 
          x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[subj_start[j]:subj_end[j]];
        
      }

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         params_to_input_p)';

      pred[subj_start[j]:subj_end[j]] = 
        x_pred[subj_start[j]:subj_end[j], 2] ./ vc_p[subj_start[j]:subj_end[j]];
        
    }

  
    for(i in 1:n_time_new){
      if(ipred[i] == 0){
        dv[i] = 0;
      }else{
        dv[i] = lognormal_rng(log(ipred[i]), sigma);
      }
    }
  }
}
