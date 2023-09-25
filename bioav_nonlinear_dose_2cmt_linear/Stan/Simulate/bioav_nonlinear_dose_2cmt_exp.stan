// First Order Absorption (oral/subcutaneous)
// Two-compartment PK Model
// Nonlinear Bioavailability - Pin the bioavailability to be 1 at
//   'reference_dose' mg. Bioavailability decreases in the form 
//   bioav = 1 - FMAX*(dose - reference_dose)^HILL/(F50^HILL + (dose - reference_dose)^HILL)
// IIV on CL, VC, Q, VP, KA (full covariance matrix)
// exponential error - DV = IPRED*exp(eps)
// Any of analytical, matrix-exponential, or general ODE solution using Torsten

functions{
  
  vector depot_2cmt_ode(real t, vector y, array[] real params, 
                        array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real q = params[3];
    real vp = params[4];
    real ka = params[5];
  
    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;

    vector[3] dydt;

    dydt[1] = -ka*y[1];                               // depot
    dydt[2] = ka*y[1] - (ke + k_cp)*y[2] + k_pc*y[3]; // central
    dydt[3] = k_cp*y[2] - k_pc*y[3];                  // peripheral
    
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
  
  real<lower = 0> TVCL;
  real<lower = 0> TVVC;
  real<lower = 0> TVQ;
  real<lower = 0> TVVP;
  real<lower = 0> TVKA;
  
  real<lower = 0, upper = 1> FMAX;
  real<lower = 0> F50;
  real<lower = 0> HILL;
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_q;
  real<lower = 0> omega_vp;
  real<lower = 0> omega_ka;
  
  corr_matrix[5] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_cl_vc, cor_cl_ka, ... and then construct the 
                     // correlation matrix in transformed data, but it's easy
                     // enough to do in R
  
  real<lower = 0> sigma;
  
  real<lower = 0> reference_dose;  // reference dose for bioavailability (lowest dose)
  array[n_subjects] real<lower = reference_dose> dose; // dose in mg for each subject
  
  int<lower = 1, upper = 3> solver; // 1 = analytical, 2 = matrix exponential, 3 = ODE
  
}
transformed data{
  
  int n_random = 5;
  int n_cmt = 3;

  vector[n_random] omega = [omega_cl, omega_vc, omega_q, omega_vp, omega_ka]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);
  
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);  // Hardcoding, but could be data or a parameter in another situation
  
  vector[n_subjects] dose_minus_reference = to_vector(dose) - reference_dose;
  
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
  vector[n_subjects] KA;
  vector[n_subjects] BIOAV;
  
  {
  
    vector[n_random] typical_values = to_vector({TVCL, TVVC, TVQ, TVVP, TVKA});
    
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
    KA = col(theta, 5);
    BIOAV = 1 - (FMAX*dose_minus_reference^HILL./(F50^HILL + dose_minus_reference^HILL));
    
    for(j in 1:n_subjects){
      
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
                           {CL[j], Q[j], VC[j], VP[j], KA[j]},
                           {BIOAV[j], 1, 1}, tlag)';
                           
      }else if(solver == 2){
        
        matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
  
        K[1, 1] = -KA[j];
        K[2, 1] = KA[j];
        K[2, 2] = -((CL[j] + Q[j])/VC[j]);
        K[2, 3] = Q[j]/VP[j];
        K[3, 2] = Q[j]/VC[j];
        K[3, 3] = -Q[j]/VP[j];
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           K, {BIOAV[j], 1, 1}, tlag)';
                           
      }else{
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(depot_2cmt_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j], Q[j], VP[j], KA[j]},
                         {BIOAV[j], 1, 1}, tlag)';
                         
      }

      ipred[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
    
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

