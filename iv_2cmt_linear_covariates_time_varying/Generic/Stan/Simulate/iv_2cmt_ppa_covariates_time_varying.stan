// IV infusion
// Two-compartment PK Model
// IIV on CL, VC, Q, and VP (full covariance matrix)
// proportional plus additive error - DV = IPRED*(1 + eps_p) + eps_a
// Any of analytical, matrix-exponential, or general ODE solution using Torsten
// Observations are generated from a normal that is truncated below at 0
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   generates values from a normal that is truncated below at 0
// Covariates - this file is generic, so all can be time-varying or constant. 
//   The key will be that the input for each covariate is length n_total and 
//   not of length n_subjects (length = n_subjects implies that covariate is 
//   constant): 
//   1) Body Weight on CL, VC, Q, VP - (wt/70)^theta
//   2) Race = Asian on VC (0/1) - exp(theta*race_asian) (obviously race is 
//        constant, but this is for an example. It'll just be constant for 
//        length n_total)
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
    
    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;
    
    vector[2] dydt;

    dydt[1] = -(ke + k_cp)*y[1] + k_pc*y[2];  // central
    dydt[2] = k_cp*y[1] - k_pc*y[2];          // peripheral
    
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
  vector<lower = 0, upper = 1>[n_total] race_asian;  // race = Asian
  vector<lower = 0>[n_total] egfr;                   // eGFR
  
  real<lower = 0> TVCL;
  real<lower = 0> TVVC;
  real<lower = 0> TVQ;
  real<lower = 0> TVVP;
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_q;
  real<lower = 0> omega_vp;
  
  real theta_cl_wt;         // allometric scaling coefficient for wt on clearance
  real theta_vc_wt;         // allometric scaling coefficient for wt on VC
  real theta_q_wt;          // allometric scaling coefficient for wt on Q
  real theta_vp_wt;         // allometric scaling coefficient for wt on VP
  real theta_vc_race_asian; // effect of race (asian) on VC 
  real theta_cl_egfr;       // effect of eGFR on clearance
  
  corr_matrix[4] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_cl_vc, cor_cl_q, ... and then construct the 
                     // correlation matrix in transformed data like is done with
                     // R_Sigma, but it's easy enough to do in R
  
  real<lower = 0> sigma_p;
  real<lower = 0> sigma_a;
  real<lower = -1, upper = 1> cor_p_a;
  
  int<lower = 1, upper = 3> solver; // 1 = analytical, 2 = matrix exponential, 3 = ODE
  
}
transformed data{
  
  int n_random = 4;
  int n_cmt = (solver == 1) ? 3 : 2;
  
  vector[n_random] omega = [omega_cl, omega_vc, omega_q, omega_vp]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);

  vector[2] sigma = [sigma_p, sigma_a]';
  matrix[2, 2] R_Sigma = rep_matrix(1, 2, 2);
  R_Sigma[1, 2] = cor_p_a;
  R_Sigma[2, 1] = cor_p_a;
  
  matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
  
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
  vector[n_total] Q;
  vector[n_total] VP;
  
  {
  
    vector[n_random] typical_values = to_vector({TVCL, TVVC, TVQ, TVVP});
    
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
      vector[n_total_subj] wt_adjustment_q = pow(wt_over_70, theta_q_wt);
      vector[n_total_subj] wt_adjustment_vp = pow(wt_over_70, theta_vp_wt);
      vector[n_total_subj] race_asian_adjustment_vc = 
                 exp(theta_vc_race_asian*race_asian[subj_start[j]:subj_end[j]]);
      vector[n_total_subj] egfr_adjustment_cl = 
                      pow(egfr[subj_start[j]:subj_end[j]] ./ 90, theta_cl_egfr);
      
      row_vector[n_random] theta_j = theta[j]; // access the parameters for subject j
      CL[subj_start[j]:subj_end[j]] = theta_j[1] .* wt_adjustment_cl .* egfr_adjustment_cl;
      VC[subj_start[j]:subj_end[j]] = theta_j[2] .* wt_adjustment_vc .* race_asian_adjustment_vc;
      Q[subj_start[j]:subj_end[j]] = theta_j[3] * wt_adjustment_q;
      VP[subj_start[j]:subj_end[j]] = theta_j[4] * wt_adjustment_vp;
      
      if(solver == 1){
        
        array[n_total_subj, n_random + 1] real params_to_input;
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
                           
      }else if(solver == 2){
        
        array[n_total] matrix[n_cmt, n_cmt] K;
        
        for(i in subj_start[j]:subj_end[j]){
          
          real ke = CL[i]/VC[i];
          real k_cp = Q[i]/VC[i];
          real k_pc = Q[i]/VP[i];
          
          K[i] = rep_matrix(0, n_cmt, n_cmt);
          
          K[i, 1, 1] = -(ke + k_cp);
          K[i, 1, 2] = k_pc;
          K[i, 2, 1] = k_cp;
          K[i, 2, 2] = -k_pc;
          
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
                           
       ipred[subj_start[j]:subj_end[j]] = 
         x_ipred[subj_start[j]:subj_end[j], 1] ./ VC[subj_start[j]:subj_end[j]];
                           
      }else{
        
        array[n_total_subj, n_random] real params_to_input;
        params_to_input[, 1] = to_array_1d(CL[subj_start[j]:subj_end[j]]);
        params_to_input[, 2] = to_array_1d(VC[subj_start[j]:subj_end[j]]);
        params_to_input[, 3] = to_array_1d(Q[subj_start[j]:subj_end[j]]);
        params_to_input[, 4] = to_array_1d(VP[subj_start[j]:subj_end[j]]);
        
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
                         params_to_input, 
                         bioav, tlag)';
                         
       ipred[subj_start[j]:subj_end[j]] = 
         x_ipred[subj_start[j]:subj_end[j], 1] ./ VC[subj_start[j]:subj_end[j]];
                         
      }

    }

    for(i in 1:n_total){
      if(ipred[i] == 0){
         dv[i] = 0;
      }else{
        real ipred_tmp = ipred[i];
        real sigma_tmp = sqrt(square(ipred_tmp) * Sigma[1, 1] + Sigma[2, 2] + 
                              2*ipred_tmp*Sigma[2, 1]);
        dv[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
        
      }
    }
  }
}

