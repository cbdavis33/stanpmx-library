// IV infusion
// Three-compartment PK Model
// IIV on CL, VC, Q1, VP1, Q2, and VP2 (full covariance matrix)
// proportional plus additive error - DV = IPRED*(1 + eps_p) + eps_a
// Either of matrix-exponential or general ODE solution using Torsten
// Observations are generated from a normal that is truncated below at 0
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   generates values from a normal that is truncated below at 0

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  vector iv_3cmt_ode(real t, vector y, array[] real params, 
                     array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real q1 = params[3];
    real vp1 = params[4];
    real q2 = params[5];
    real vp2 = params[6];
    
    real ke = cl/vc;
    real k_cp1 = q1/vc;
    real k_p1c = q1/vp1;
    real k_cp2 = q2/vc;
    real k_p2c = q2/vp2;
    
    vector[3] dydt;

    dydt[1] = -(ke + k_cp1 + k_cp2)*y[1] + k_p1c*y[2] + k_p2c*y[3];  // central
    dydt[2] = k_cp1*y[1] - k_p1c*y[2];                               // peripheral 1
    dydt[3] = k_cp2*y[1] - k_p2c*y[3];                               // peripheral 2
    
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
  real<lower = 0> TVQ1;
  real<lower = 0> TVVP1;
  real<lower = 0> TVQ2;
  real<lower = 0> TVVP2;
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_q1;
  real<lower = 0> omega_vp1;
  real<lower = 0> omega_q2;
  real<lower = 0> omega_vp2;
  
  corr_matrix[6] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_cl_vc, cor_cl_q1, ... and then construct the 
                     // correlation matrix in transformed data like is done with
                     // R_Sigma, but it's easy enough to do in R
  
  real<lower = 0> sigma_p;
  real<lower = 0> sigma_a;
  real<lower = -1, upper = 1> cor_p_a;
  
  int<lower = 1, upper = 2> solver; // 1 = matrix exponential, 2 = rk45
  
}
transformed data{
  
  int n_random = 6;
  int n_cmt = 3;
  
  vector[n_random] omega = [omega_cl, omega_vc, omega_q1, omega_vp1, 
                                                omega_q2, omega_vp2]';
  
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
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q1;
  vector[n_subjects] VP1;
  vector[n_subjects] Q2;
  vector[n_subjects] VP2;
  
  {
  
    vector[n_random] typical_values = to_vector({TVCL, TVVC, TVQ1, TVVP1,
                                                             TVQ2, TVVP2});
    
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
    Q1 = col(theta, 3);
    VP1 = col(theta, 4);
    Q2 = col(theta, 5);
    VP2 = col(theta, 6);
    
    for(j in 1:n_subjects){
                           
      if(solver == 1){
        
        real ke = CL[j]/VC[j];
        real k_cp1 = Q1[j]/VC[j];
        real k_p1c = Q1[j]/VP1[j];
        real k_cp2 = Q2[j]/VC[j];
        real k_p2c = Q2[j]/VP2[j];
        
        matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
        K[1, 1] = -(ke + k_cp1 + k_cp2);
        K[1, 2] = k_p1c;
        K[1, 3] = k_p2c;
        K[2, 1] = k_cp1;
        K[2, 2] = -k_p1c;
        K[3, 1] = k_cp2;
        K[3, 3] = -k_p2c;
        
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
          pmx_solve_rk45(iv_3cmt_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j], Q1[j], VP1[j], Q2[j], VP2[j]}, 
                         bioav, tlag)';
                         
      }

      ipred[subj_start[j]:subj_end[j]] = 
                      x_ipred[subj_start[j]:subj_end[j], 1] ./ VC[j];
    
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

