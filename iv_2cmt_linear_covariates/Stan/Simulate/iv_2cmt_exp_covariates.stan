// IV infusion
// One-compartment PK Model
// IIV on CL, VC, Q, and VP (full covariance matrix)
// exponential error - DV = IPRED*exp(eps)
// Any of analytical, matrix-exponential, or general ODE solution using Torsten
// Covariates: 
//   1) Body Weight on CL, VC, Q, and VP - (wt/70)^theta
//   2) Sex on VC (0 = male, 1 = female) - exp(theta*sex)
//   3) eGFR on CL (continuous) - (eGFR/90)^theta
//   4) Race on VC - exp(theta_vc_{race}*I(race == {race}))

functions{
  
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
  
  array[n_subjects] real<lower = 0> wt;                    // baseline bodyweight (kg)
  array[n_subjects] int<lower = 0, upper = 1> sex;         // sex
  array[n_subjects] real<lower = 0> egfr;                  // eGFR
  int<lower = 2> n_races;                                  // number of unique races
  array[n_subjects] int<lower = 1, upper = n_races> race;  // race

  
  real<lower = 0> TVCL;
  real<lower = 0> TVVC;
  real<lower = 0> TVQ;
  real<lower = 0> TVVP;
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_q;
  real<lower = 0> omega_vp;
  
  real theta_cl_wt;    // allometric scaling coefficient for wt on clearance
  real theta_vc_wt;    // allometric scaling coefficient for wt on VC
  real theta_q_wt;     // allometric scaling coefficient for wt on Q
  real theta_vp_wt;    // allometric scaling coefficient for wt on VP
  real theta_vc_sex;   // effect of sex on VC 
  real theta_cl_egfr;  // effect of eGFR on clearance
  real theta_vc_race2; // effect of race == 2 relative to race == 1 on VC
  real theta_vc_race3; // effect of race == 3 relative to race == 1 on VC
  real theta_vc_race4; // effect of race == 4 relative to race == 1 on VC
  
  corr_matrix[4] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_cl_vc and then construct the 
                     // correlation matrix in transformed data like is done with
                     // R_Sigma, but it's easy enough to do in R
  
  real<lower = 0> sigma;
  
  int<lower = 1, upper = 3> solver; // 1 = analytical, 2 = matrix exponential, 3 = ODE
  
}
transformed data{
  
  int n_random = 4;
  int n_cmt = (solver == 1) ? 3 : 2;
  
  vector[n_races] theta_vc_race = [0, theta_vc_race2, theta_vc_race3, 
                                   theta_vc_race4]'; 
  
  vector[n_random] omega = [omega_cl, omega_vc, omega_q, omega_vp]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);  // Hardcoding, but could be data or a parameter in another situation
  
}
model{
  
}
generated quantities{
  
  vector[n_total] ipred; // concentration with no residual error
  vector[n_total] dv;    // concentration with residual error
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;
  
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

    CL = col(theta, 1);
    VC = col(theta, 2);
    Q = col(theta, 3);
    VP = col(theta, 4);
    
    for(j in 1:n_subjects){
      
      real wt_over_70 = wt[j]/70;
      real wt_adjustment_cl = wt_over_70^theta_cl_wt;
      real wt_adjustment_vc = wt_over_70^theta_vc_wt;
      real wt_adjustment_q = wt_over_70^theta_q_wt;
      real wt_adjustment_vp = wt_over_70^theta_vp_wt;
      real sex_adjustment_vc = exp(theta_vc_sex*sex[j]);
      real egfr_adjustment_cl = (egfr[j]/90)^theta_cl_egfr;
      real race_adjustment_vc = exp(theta_vc_race[race[j]]);
      
      row_vector[n_random] theta_j = theta[j]; // access the parameters for subject j
      CL[j] = theta_j[1] * wt_adjustment_cl * egfr_adjustment_cl;
      VC[j] = theta_j[2] * wt_adjustment_vc * race_adjustment_vc * sex_adjustment_vc;
      Q[j] = theta_j[3] * wt_adjustment_q;
      VP[j] = theta_j[4] * wt_adjustment_vp;
      
      if(solver == 1){
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           {CL[j], Q[j], VC[j], VP[j], 0})';
                           
      }else if(solver == 2){
        
        real ke = CL[j]/VC[j];
        real k_cp = Q[j]/VC[j];
        real k_pc = Q[j]/VP[j];
        
        matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
        K[1, 1] = -(ke + k_cp);
        K[1, 2] = k_pc;
        K[2, 1] = k_cp;
        K[2, 2] = -k_pc;
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           K, bioav, tlag)';
                           
      }else{
        
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
                         bioav, tlag)';
                         
      }

      ipred[subj_start[j]:subj_end[j]] = 
                      x_ipred[subj_start[j]:subj_end[j], (n_cmt - 1)] ./ VC[j];
    
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


