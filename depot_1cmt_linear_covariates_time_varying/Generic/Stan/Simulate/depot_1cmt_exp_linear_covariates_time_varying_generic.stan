// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// IIV on CL, VC, KA (full covariance matrix)
// exponential error - DV = IPRED*exp(eps)
// Any of analytical, matrix-exponential, or general ODE solution using Torsten
// Covariates - this file is generic, so all can be time-varying or constant. 
//   The key will be that the input for each covariate is length n_total and 
//   not of length n_subjects (length = n_subjects implies that covariate is 
//   constant): 
//   1) Body Weight on CL and VC - (wt/70)^theta
//   2) Concomitant administration of protein pump inhibitors (CMPPI) 
//      on KA (0/1) - exp(theta*cmppi)
//   3) eGFR on CL (continuous) - (eGFR/90)^theta

functions{
  
  vector depot_1cmt_ode(real t, vector y, array[] real params, 
                        array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real ke = cl/vc;
    
    vector[2] dydt;

    dydt[1] = -ka*y[1];           // depot
    dydt[2] = ka*y[1] - ke*y[2];  // central
    
    return dydt;
  }
  
}
data{
  
  int n_subjects;
  int n_total;                  
  
  array[n_total] real amt;
  array[n_total] int cmt;
  array[n_total] int evid;
  array[n_total] real rate;
  array[n_total] real ii;
  array[n_total] int addl;
  array[n_total] int ss;
  array[n_total] real time;
  
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;
  
  vector<lower = 0>[n_total] wt;                     // bodyweight (kg)
  vector<lower = 0, upper = 1>[n_total] cmppi;       // cmppi
  vector<lower = 0>[n_total] egfr;                   // eGFR
  
  real<lower = 0> TVCL;
  real<lower = 0> TVVC;
  real<lower = 0> TVKA;
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_ka;
  
  real theta_cl_wt;    // allometric scaling coefficient for wt on clearance
  real theta_vc_wt;    // allometric scaling coefficient for wt on VC
  real theta_ka_cmppi; // effect of CMPPI on KA 
  real theta_cl_egfr;  // effect of eGFR on clearance
  
  corr_matrix[3] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_cl_vc, cor_cl_ka, ... and then construct the 
                     // correlation matrix in transformed data, but it's easy
                     // enough to do in R
  
  real<lower = 0> sigma;
  
  int<lower = 1, upper = 3> solver; // 1 = analytical, 2 = matrix exponential, 3 = ODE
  
}
transformed data{
  
  int n_random = 3;
  int n_cmt = 2;

  vector[n_random] omega = [omega_cl, omega_vc, omega_ka]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);

  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);  // Hardcoding, but could be data or a parameter in another situation
  
}
model{
  
}
generated quantities{
  
  vector[n_total] ipred; // concentration with no residual error
  vector[n_total] dv;    // concentration with residual error
  
  vector[n_total] CL;
  vector[n_total] VC;
  vector[n_total] KA;
  
  {
  
    vector[n_random] typical_values = to_vector({TVCL, TVVC, TVKA});
    
    matrix[n_random, n_subjects] eta;   
    matrix[n_subjects, n_random] theta; 
  
    matrix[n_total, n_cmt] x_ipred;
    
    for(i in 1:n_subjects){
      eta[, i] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                           diag_pre_multiply(omega, L));
    }
    theta = (rep_matrix(typical_values, n_subjects) .* exp(eta))';
    
    for(j in 1:n_subjects){
      
      int n_total_subj = subj_end[j] - subj_start[j] + 1;
      vector[n_total_subj] wt_over_70 = wt[subj_start[j]:subj_end[j]] ./ 70;
      vector[n_total_subj] wt_adjustment_cl = pow(wt_over_70, theta_cl_wt);
      vector[n_total_subj] wt_adjustment_vc = pow(wt_over_70, theta_vc_wt);
      vector[n_total_subj] cmppi_adjustment_ka = 
                      exp(theta_ka_cmppi*cmppi[subj_start[j]:subj_end[j]]);
      vector[n_total_subj] egfr_adjustment_cl = 
                      pow(egfr[subj_start[j]:subj_end[j]] ./ 90, theta_cl_egfr);
      
      row_vector[n_random] theta_j = theta[j]; // access the parameters for subject j
      CL[subj_start[j]:subj_end[j]] = theta_j[1] .* wt_adjustment_cl .* egfr_adjustment_cl;
      VC[subj_start[j]:subj_end[j]] = theta_j[2] * wt_adjustment_vc;
      KA[subj_start[j]:subj_end[j]] = theta_j[3] * cmppi_adjustment_ka;
 
      if(solver == 1){
        
        array[n_total_subj, n_random] real params_to_input;
        params_to_input[, 1] = to_array_1d(CL[subj_start[j]:subj_end[j]]);
        params_to_input[, 2] = to_array_1d(VC[subj_start[j]:subj_end[j]]);
        params_to_input[, 3] = to_array_1d(KA[subj_start[j]:subj_end[j]]);
        
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
                           
      }else if(solver == 2){
        
        array[n_total] matrix[n_cmt, n_cmt] K;
        
        for(i in subj_start[j]:subj_end[j]){
          
          real ke = CL[i]/VC[i];
          
          K[i] = rep_matrix(0, n_cmt, n_cmt);
          
          K[i, 1, 1] = -KA[i];
          K[i, 2, 1] = KA[i];
          K[i, 2, 2] = -ke;
          
        }
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           K[subj_start[j]:subj_end[j]], 
                           bioav, tlag)';
                           
      }else{
        
        array[n_total_subj, n_random] real params_to_input;
        params_to_input[, 1] = to_array_1d(CL[subj_start[j]:subj_end[j]]);
        params_to_input[, 2] = to_array_1d(VC[subj_start[j]:subj_end[j]]);
        params_to_input[, 3] = to_array_1d(KA[subj_start[j]:subj_end[j]]);
        
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
                         params_to_input, 
                         bioav, tlag)';
                         
      }

      ipred[subj_start[j]:subj_end[j]] = 
         x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[subj_start[j]:subj_end[j]];
    
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
