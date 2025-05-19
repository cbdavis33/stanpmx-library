// IV infusion
// Two-compartment PK Model
// IIV on CL, VC, Q, and VP (full covariance matrix)
// exponential error - DV = IPRED*exp(eps)
// Output includes individual half-life, AUC since 0 for every timepoint, and 
//   AUC between t1 and t2 (like a dosing interval). If you want Cmax and Tmax, 
//   make sure to simulate the end of the infusion and pull that time and 
//   concentration out.

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

    dydt[1] = -(ke + k_cp)*y[1] + k_pc*y[2];   // central
    dydt[2] = k_cp*y[1] - k_pc*y[2];           // peripheral
    dydt[3] = y[1];                            // AUC
    dydt[4] = t >= t_1 && t <= t_2 ? y[1] : 0; // AUC_t_1-t_2
    
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
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_q;
  real<lower = 0> omega_vp;
  
  corr_matrix[4] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_cl_vc and then construct the 
                     // correlation matrix in transformed data like is done with
                     // R_Sigma, but it's easy enough to do in R
  
  real<lower = 0> sigma;
  
  real<lower = 0> t_1;
  real<lower = t_1> t_2;
  
}
transformed data{
  
  int n_random = 4;
  int n_cmt = 4;
  
  vector[n_random] omega = [omega_cl, omega_vc, omega_q, omega_vp]';
  
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
  vector[n_subjects] t_half_alpha;    // alpha half-life
  vector[n_subjects] t_half_terminal; // terminal half-life
  
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
    
    vector[n_subjects] alpha = 0.5*(CL./VC + Q./VC + Q./VP + 
                            sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
    vector[n_subjects] beta = 0.5*(CL./VC + Q./VC + Q./VP - 
                            sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
    
    t_half_alpha = log(2)/alpha;
    t_half_terminal = log(2)/beta;
    
    for(j in 1:n_subjects){

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
                                
      auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 3] ./ VC[j];
      
      auc_t1_t2[j] = max(x_ipred[subj_start[j]:subj_end[j], 4]) / VC[j];
    
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



