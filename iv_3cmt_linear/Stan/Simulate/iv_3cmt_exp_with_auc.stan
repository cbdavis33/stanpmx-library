// IV infusion
// Three-compartment PK Model
// IIV on CL, VC, Q1, VP1, Q2, and VP2 (full covariance matrix)
// proportional plus additive error - DV = IPRED*(1 + eps_p) + eps_a
// General ODE solution using Torsten
// Output includes individual AUC since 0 for every timepoint, AUC between 
//   t1 and t2 (like a dosing interval), and each elimination half-life. If you 
//   want Cmax and Tmax, make sure to simulate the end of the infusion and pull 
//   that time and concentration out.

functions{
  
  vector iv_3cmt_ode(real t, vector y, array[] real params, 
                     array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real q1 = params[3];
    real vp1 = params[4];
    real q2 = params[5];
    real vp2 = params[6];
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real ke = cl/vc;
    real k_cp1 = q1/vc;
    real k_p1c = q1/vp1;
    real k_cp2 = q2/vc;
    real k_p2c = q2/vp2;
    
    vector[5] dydt;

    dydt[1] = -(ke + k_cp1 + k_cp2)*y[1] + k_p1c*y[2] + k_p2c*y[3];  // central
    dydt[2] = k_cp1*y[1] - k_p1c*y[2];                               // peripheral 1
    dydt[3] = k_cp2*y[1] - k_p2c*y[3];                               // peripheral 2
    dydt[4] = y[1];                                                  // AUC
    dydt[5] = t >= t_1 && t <= t_2 ? y[1] : 0;                       // AUC_t_1-t_2
    
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
  
  real<lower = 0> sigma;
  
  real<lower = 0> t_1;
  real<lower = t_1> t_2;
  
}
transformed data{
  
  int n_random = 6;
  int n_cmt = 5;
  
  vector[n_random] omega = [omega_cl, omega_vc, omega_q1, omega_vp1, 
                                                omega_q2, omega_vp2]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);

  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);  // Hardcoding, but could be data or a parameter in another situation
  
  array[1, 2] real x_r = {{t_1, t_2}};
  
}
model{
  
}
generated quantities{
  
  vector[n_total] ipred; // concentration with no residual error
  vector[n_total] dv;    // concentration with residual error
  
  vector[n_total] auc;                // AUC 
  vector[n_subjects] auc_t1_t2;       // AUC from t1 up to t2
  vector[n_subjects] t_half_alpha;    // half_life of the first elimination phase
  vector[n_subjects] t_half_beta;     // half_life of the second elimination phase
  vector[n_subjects] t_half_terminal; // half_life of the final elimination phase
  
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
        
      real ke = CL[j]/VC[j];
      real k_cp1 = Q1[j]/VC[j];
      real k_p1c = Q1[j]/VP1[j];
      real k_cp2 = Q2[j]/VC[j];
      real k_p2c = Q2[j]/VP2[j];
      
      real jay = k_cp1 + ke + k_p1c + k_p2c + k_cp2;
      real kay = k_cp1*k_p2c + ke*k_p1c + ke*k_p2c + k_p1c*k_p2c + k_cp2*k_p1c;
      real ell = ke*k_p1c*k_p2c;
      real m = (3*kay - square(jay))/3;
      real n = (2*jay^3 - 9*jay*kay + 27*ell)/27;
      real Q = square(n)/4 + m^3/27;
      real alpha = sqrt(-Q);
      real beta = -n/2;
      real rho = hypot(beta, alpha);  // sqrt(square(beta) + square(alpha))
      real delta = atan2(alpha, beta);
      
      real lambda_1 = jay/3 + cbrt(rho)*(cos(delta/3) + sqrt(3)*sin(delta/3));
      real lambda_2 = jay/3 + cbrt(rho)*(cos(delta/3) - sqrt(3)*sin(delta/3));
      real lambda_3 = jay/3 - 2*cbrt(rho)*cos(delta/3);
        
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
                       bioav, tlag, x_r)';
                         
      ipred[subj_start[j]:subj_end[j]] = 
                      x_ipred[subj_start[j]:subj_end[j], 1] ./ VC[j];
                      
      auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 4] ./ VC[j];
      
      auc_t1_t2[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]) / VC[j];
      
      t_half_alpha[j] = log(2)/lambda_1;
      t_half_beta[j] = log(2)/lambda_2;
      t_half_terminal[j] = log(2)/lambda_3;
    
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
