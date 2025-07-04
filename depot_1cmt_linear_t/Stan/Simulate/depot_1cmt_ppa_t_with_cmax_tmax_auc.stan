// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// IIV on CL, VC, KA (full covariance matrix)
// proportional plus additive error with a Student's-t distribution - 
//   DV = IPRED*(1 + eps_p) + eps_a
// User's choice of analytical, matrix-exponential, ODE, or ODE with AUC, Cmax, 
//   Tmax
// Observations are generated from a Student's-t that is truncated below at 0
// Since we have a Student's-t distribution on the error, but the DV must 
//   be > 0, it generates values from a t that is truncated below at 0
// If solver == 4 (ODE with AUC, Cmax, Tmax), output includes individual Cmax 
//   over the whole time period, Tmax between t1 and t2, AUC since 0 for every 
//   timepoint, and AUC between t1 and t2 (like a dosing interval)

functions{
  
  vector n_t_q_inner_coeffs(int n, real nu) {
    if (n <= 0)
      reject("n must be positive; found n = ", n);
  
  real g = sqrt(nu/2) * exp(lgamma(nu/2) - lgamma((nu + 1)/2));
    vector[3] x0 = [
      g,
      ((nu + 1)*g^3 - nu*g)/(6*nu),
      ((7*nu^2 + 8*nu + 1)*g^5 - 10*(nu^2 + nu)*g^3 + 3*nu^2*g)/(120*nu^2)
    ]';
  
    if(n < 4)
      return x0[1:n];
  
    vector[n] x = append_row(x0, rep_vector(0, n - 3));
  
    for(i in 2:(n - 2)) {
      x[i + 2] = -(2*i + 1)*x[i + 1];
      
      for(l in 0:i) {
        for(m in 0:(i - l)) {
          x[i + 2] = x[i + 2] + ((1 + 1/nu)*(2*l + 1)*(2*m + 1) - 
                       2*m*(2*m + 1)/nu)*x[i - l - m + 1]*x[l + 1]*x[m + 1];
        }
      }
      for (l in 0:(i - 1)) {
        for (m in 0:(i - l - 1)) {
          x[i + 2] = x[i + 2] - (1/nu)*(2*m + 1)*x[i - l - m]*x[l + 1]*x[m + 1];
        }
      }
      
      x[i + 2] = x[i + 2]/(2*i + 3)/(2*i + 2);
    }
  
    return x;
}

  real n_t_q(real z, real nu) {
    if(abs(z) < 2*sqrt(nu)) {
      vector[11] a = n_t_q_inner_coeffs(11, nu);
      for (i in 0:10)
        a[i + 1] *= z^(2*i + 1);
      return sum(a);
    } else {
      real w = (1 - std_normal_cdf(abs(z)))*nu*sqrt(pi())*exp(lgamma(nu/2) - lgamma((nu + 1)/2));
      real sign_z = z > 0 ? 1 : z < 0 ? -1 : 0;
      return sign_z*sqrt(nu)*w^(-1/nu)*(1 - (nu + 1)/(2*(nu + 2))*w^(2/nu));
    }
  }

  real student_t_alternate_rng(real nu, real mu, real sigma) {
    return mu + sigma * n_t_q(std_normal_rng(), nu);
  }

  real student_t_lub_rng(real nu, real mu, real sigma, real lb, real ub) {
    real p_lb = student_t_cdf(lb | nu, mu, sigma);
    real p_ub = student_t_cdf(ub | nu, mu, sigma);
    real u = uniform_rng(p_lb, p_ub);
    real z = inv_Phi(u);
    return mu + sigma * n_t_q(z, nu);
  }
  
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
  
  vector depot_1cmt_with_auc_ode(real t, vector y, array[] real params, 
                                 array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real ke = cl/vc;
    real slope = ka*y[1] - ke*y[2];
    real x = slope > 0 && y[2]/vc > y[5] ? slope/vc : 0;
    real z = t <= t_1 || (slope > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    vector[6] dydt;

    dydt[1] = -ka*y[1];                            // depot
    dydt[2] = slope;                               // central
    dydt[3] = y[2];                                // AUC
    dydt[4] = t >= t_1 && t <= t_2 ? y[2] : 0;     // AUC_t_1-t_2
    dydt[5] = x;                                   // C_max
    dydt[6] = z;                                   // t_max
    
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
  real<lower = 0> TVKA;
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_ka;
  
  corr_matrix[3] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_cl_vc, cor_cl_ka, ... and then construct the 
                     // correlation matrix in transformed data, but it's easy
                     // enough to do in R
  
  real<lower = 0> sigma_p;
  real<lower = 0> sigma_a;
  real<lower = -1, upper = 1> cor_p_a;
  
  real<lower = 2> nu;
  
  int<lower = 1, upper = 4> solver; // 1 = analytical, 2 = mat-exp, 3 = rk45, 4 = rk45 with AUC, Cmax, Tmax
  
  real<lower = 0> t_1;
  real<lower = t_1> t_2;
  
}
transformed data{
  
  int n_random = 3;
  int n_cmt = (solver == 4) ? 6 : 2;

  vector[n_random] omega = [omega_cl, omega_vc, omega_ka]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);

  vector[2] sigma = [sigma_p, sigma_a]';
  matrix[2, 2] R_Sigma = rep_matrix(1, 2, 2);
  R_Sigma[1, 2] = cor_p_a;
  R_Sigma[2, 1] = cor_p_a;
  
  matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);  // Hardcoding, but could be data or a parameter in another situation
  
  array[1, 2] real x_r = {{t_1, t_2}};
  
}
model{
  
}
generated quantities{
  
  vector[n_total] ipred; // concentration with no residual error
  vector[n_total] dv;    // concentration with residual error
  
  vector[n_total] auc;            // AUC 
  vector[n_subjects] auc_t1_t2;   // AUC from t1 up to t2
  vector[n_subjects] c_max;       // Cmax
  vector[n_subjects] t_max;       // Tmax
  vector[n_subjects] t_half;      // half-life
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  
  {
  
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA});
    
    matrix[n_subjects, n_random] eta;   
    matrix[n_subjects, n_random] theta; 
  
    matrix[n_total, n_cmt] x_ipred;
    
    for(i in 1:n_subjects){
      eta[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                           diag_pre_multiply(omega, L))';
    }
    theta = (rep_matrix(typical_values, n_subjects) .* exp(eta));

    CL = col(theta, 1);
    VC = col(theta, 2);
    KA = col(theta, 3);

    for(j in 1:n_subjects){
      
      t_half[j] = log(2)/(CL[j]/VC[j]);
      
      if(solver == 1){
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           {CL[j], VC[j], KA[j]},
                           bioav, tlag)';
        
      }else if(solver == 2){
        
        matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
  
        K[1, 1] = -KA[j];
        K[2, 1] = KA[j];
        K[2, 2] = -CL[j]/VC[j];
        
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
        
      }else if(solver == 3){
        
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
                         {CL[j], VC[j], KA[j]}, 
                         bioav, tlag, x_r)';
        
      }else{
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(depot_1cmt_with_auc_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j], KA[j]}, 
                         bioav, tlag, x_r)';
        
        auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 3] ./ VC[j];
      
        auc_t1_t2[j] = max(x_ipred[subj_start[j]:subj_end[j], 4]) / VC[j];
        c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]);
        t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 6]) - t_1;
        
      }
      
      ipred[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
    
    }

    for(i in 1:n_total){
      if(ipred[i] == 0){
         dv[i] = 0;
      }else{
        real ipred_tmp = ipred[i];
        real sigma_tmp = sqrt(square(ipred_tmp) * Sigma[1, 1] + Sigma[2, 2] + 
                              2*ipred_tmp*Sigma[2, 1]);
        dv[i] = student_t_lub_rng(nu, ipred_tmp, sigma_tmp, 0.0, 
                                  positive_infinity());
        
      }
    }
  }
}

