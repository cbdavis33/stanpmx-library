// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// IIV on CL, VC, and Ka (full covariance matrix)
// proportional plus additive error - DV = IPRED*(1 + eps_p) + eps_a
// General ODE solution using Torsten to get out individual estimates of AUC, 
//   Cmax, Tmax, ... is an option, but the user can choose to only calculate 
//   IPRED and DV
// Predictions are generated from a normal that is truncated below at 0
// Covariates - this file is generic, so all can be time-varying or constant. 
//   The key will be that the input for each covariate is length n_total and 
//   not of length n_subjects (length = n_subjects implies that covariate is 
//   constant): 
//   1) Body Weight on CL and VC - (wt/70)^theta
//   2) Concomitant administration of protein pump inhibitors (CMPPI) 
//      on KA (0/1) - exp(theta*cmppi)
//   3) eGFR on CL (continuous) - (eGFR/90)^theta

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  vector depot_1cmt_ode(real t, vector y, array[] real params, 
                        array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real ke = cl/vc;
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
   
    real slope = ka*y[1] - ke*y[2];
    real x = slope > 0 && y[2]/vc > y[4] ? slope/vc : 0;
    real z = t <= t_1 || (slope > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    vector[5] dydt;

    dydt[1] = -ka*y[1];                            // depot
    dydt[2] = slope;                               // central
    dydt[3] = t >= t_1 && t <= t_2 ? y[2]/vc : 0;  // AUC_t_1-t_2
    dydt[4] = x;                                   // C_max
    dydt[5] = z;                                   // t_max
    
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
  
  vector<lower = 0>[n_time_new] wt;                     // bodyweight (kg)
  vector<lower = 0, upper = 1>[n_time_new] cmppi;       // cmppi
  vector<lower = 0>[n_time_new] egfr;                   // eGFR
  
  real<lower = 0> t_1;   // Time at which to start SS calculations (AUC_ss, C_max_ss, ...)
  real<lower = t_1> t_2; // Time at which to end SS calculations (AUC_ss, C_max_ss, ...)
  
  int <lower = 0, upper = 0> want_auc_cmax; // 0 => only calculate concentrations
                                            // 1 => concentrations and auc, c_max, ...
                                            // For now, the c_max, t_max, and auc 
                                            // calculations aren't working as expected
                                            // with the time-varying parameters,
                                            // so this must be 0
 
}
transformed data{ 
  
  int n_random = 3; // Number of random effects
  int n_cmt = want_auc_cmax ? 5 : 2; // Number of compartments - depot, central (AUC_ss, Cmax_ss, Tmax_ss))
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt);
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  array[1, 2] real x_r = {{t_1, t_2}};

}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = TVCL/TVVC> TVKA;
  
  real theta_cl_wt;
  real theta_vc_wt;
  real theta_ka_cmppi;
  real theta_cl_egfr;
  
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
  vector[want_auc_cmax ? n_subjects : 0] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[want_auc_cmax ? n_subjects : 0] c_max;   // Cmax between t1 and t2 (c_max_ss)
  vector[want_auc_cmax ? n_subjects : 0] t_max;   // Tmax between t1 and t2, then subtract off t1
  
  vector[n_time_new] CL;
  vector[n_time_new] VC;
  vector[n_time_new] KA;
  vector[n_time_new] KE;
 
  {
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA});
    
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);

    matrix[n_subjects_new, n_random] eta_new;
    matrix[n_subjects_new, n_random] theta_new;
    matrix[n_time_new, n_cmt] x_ipred;
    matrix[n_time_new, 2] x_pred;
    
    real sigma_sq_p = square(sigma[1]);
    real sigma_sq_a = square(sigma[2]);
    
    real sigma_p_a = Sigma[1, 2];

    for(i in 1:n_subjects_new){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects_new) .* exp(eta_new));
    
    vector[n_time_new] cl_p;
    vector[n_time_new] vc_p;
    vector[n_time_new] ka_p;

    for(j in 1:n_subjects_new){
      
      int n_total_subj = subj_end[j] - subj_start[j] + 1;
      vector[n_total_subj] wt_over_70 = wt[subj_start[j]:subj_end[j]] ./ 70;
      vector[n_total_subj] wt_adjustment_cl = pow(wt_over_70, theta_cl_wt);
      vector[n_total_subj] wt_adjustment_vc = pow(wt_over_70, theta_vc_wt);
      vector[n_total_subj] cmppi_adjustment_ka = 
                      exp(theta_ka_cmppi*cmppi[subj_start[j]:subj_end[j]]);
      vector[n_total_subj] egfr_adjustment_cl = 
                      pow(egfr[subj_start[j]:subj_end[j]] ./ 90, theta_cl_egfr);
                      
      array[n_total_subj, n_random] real params_to_input;
      array[n_total_subj, n_random] real params_to_input_p;
                      
      row_vector[n_random] theta_j = theta_new[j]; // access the parameters for subject j
      CL[subj_start[j]:subj_end[j]] = theta_j[1] .* wt_adjustment_cl .* egfr_adjustment_cl;
      VC[subj_start[j]:subj_end[j]] = theta_j[2] * wt_adjustment_vc;
      KA[subj_start[j]:subj_end[j]] = theta_j[3] * cmppi_adjustment_ka;
      
      cl_p[subj_start[j]:subj_end[j]] = TVCL .* wt_adjustment_cl .* egfr_adjustment_cl;
      vc_p[subj_start[j]:subj_end[j]] = TVVC * wt_adjustment_vc;
      ka_p[subj_start[j]:subj_end[j]] = TVKA * cmppi_adjustment_ka;
      
      params_to_input[, 1] = to_array_1d(CL[subj_start[j]:subj_end[j]]);
      params_to_input[, 2] = to_array_1d(VC[subj_start[j]:subj_end[j]]);
      params_to_input[, 3] = to_array_1d(KA[subj_start[j]:subj_end[j]]);
      
      params_to_input_p[, 1] = to_array_1d(cl_p[subj_start[j]:subj_end[j]]);
      params_to_input_p[, 2] = to_array_1d(vc_p[subj_start[j]:subj_end[j]]);
      params_to_input_p[, 3] = to_array_1d(ka_p[subj_start[j]:subj_end[j]]);
      
      if(want_auc_cmax == 1){
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(depot_1cmt_ode,
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
                         
        auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], 3]);
        c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 4]);
        t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]) - t_1;
                         
      }else{
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           params_to_input)';
        
      }
                      
      ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[subj_start[j]:subj_end[j]];

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
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
        real ipred_tmp = ipred[i];
        real sigma_tmp = sqrt(square(ipred_tmp) * sigma_sq_p + sigma_sq_a + 
                              2*ipred_tmp*sigma_p_a);
        dv[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
      }
    }
  }
}
