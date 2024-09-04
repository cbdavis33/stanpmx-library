// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model 
// Inhibitory Effect Compartment PD Model
// IIV on CL, VC, and KA (full covariance matrix)
// IIV on E0, KE0, and EC50 (full covariance matrix) for PD. EMAX and HILL are 
//   fixed to be 1
// For PD - IPRED = E0*(1 - emax*conc_effcmt^hill/(ec50^hill + conc_effcmt^hill))
// proportional plus additive error for PK - DV = IPRED*(1 + eps_p) + eps_a
// proportional plus additive error for PD - DV = IPRED*(1 + eps_p_pd) + eps_a_pd
// Users choice of matrix-exponential, general ODE, or coupled ODE solution 
//   using Torsten
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
  
  vector depot_1cmt_effcmt_ode(real t, vector y, array[] real params, 
                               array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real ke0 = params[4];
    
    real ke = cl/vc;
    
    vector[3] dydt;

    dydt[1] = -ka*y[1];
    dydt[2] = ka*y[1] - ke*y[2];
    dydt[3] = ke0*(y[2] - y[3]);
    
    return dydt;
  }
  
  vector depot_1cmt_effcmt_ode_coupled(real t, vector y, vector y_pk, 
                                       array[] real params, array[] real x_r, 
                                       array[] int x_i){
    
    real vc = params[2];

    real ke0 = params[4];
    
    vector[1] dydt;

    dydt[1] = ke0*(y_pk[2] - y[1]);
    
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
  
  real<lower = 0> TVE0;
  real<lower = 0> TVKE0;
  real<lower = 0> TVEC50;
  
  real<lower = 0> omega_e0;
  real<lower = 0> omega_ke0;
  real<lower = 0> omega_ec50;
  
  corr_matrix[3] R;     // Correlation matrix before transforming to Omega.
                        // Can in theory change this to having inputs for
                        // cor_cl_vc, cor_cl_ka, ... and then construct the 
                        // correlation matrix in transformed data, but it's easy
                        // enough to do in R
                        
  corr_matrix[3] R_pd;  // Correlation matrix before transforming to Omega.
                        // Can in theory change this to having inputs for
                        // cor_e0_ke0, cor_ke0_ec50, ... and then construct the 
                        // correlation matrix in transformed data, but it's easy
                        // enough to do in R
  
  real<lower = 0> sigma_p;
  real<lower = 0> sigma_a;
  real<lower = -1, upper = 1> cor_p_a;
  
  real<lower = 0> sigma_p_pd;
  real<lower = 0> sigma_a_pd;
  real<lower = -1, upper = 1> cor_p_a_pd;
  
  int<lower = 1, upper = 3> solver; // 1 = mat-exp, 2 = General ODE, 3 = Coupled ODE 
  
}
transformed data{
  
  int n_random = 3;
  int n_random_pd = 3;
  
  int n_cmt = 2;     // number of ODEs in PK model (depot, central)
  int n_cmt_pd = 1;  // number of ODEs in PD system

  vector[n_random] omega = [omega_cl, omega_vc, omega_ka]';
                               
  vector[n_random_pd] omega_pd = [omega_e0, omega_ke0, omega_ec50]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);
  
  vector[2] sigma = [sigma_p, sigma_a]';
  matrix[2, 2] R_Sigma = rep_matrix(1, 2, 2);
  R_Sigma[1, 2] = cor_p_a;
  R_Sigma[2, 1] = cor_p_a;
  
  matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
  
  matrix[n_random_pd, n_random_pd] L_pd = cholesky_decompose(R_pd);
  
  vector[2] sigma_pd = [sigma_p_pd, sigma_a_pd]';
  matrix[2, 2] R_Sigma_pd = rep_matrix(1, 2, 2);
  R_Sigma_pd[1, 2] = cor_p_a_pd;
  R_Sigma_pd[2, 1] = cor_p_a_pd;
  
  matrix[2, 2] Sigma_pd = quad_form_diag(R_Sigma_pd, sigma_pd);
  
  array[n_cmt + n_cmt_pd] real bioav = rep_array(1.0, n_cmt + n_cmt_pd); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt + n_cmt_pd] real tlag = rep_array(0.0, n_cmt + n_cmt_pd);  // Hardcoding, but could be data or a parameter in another situation
  
  real emax = 1.0;
  real hill = 1.0;
  
}
model{
  
}
generated quantities{
  
  vector[n_total] ipred; // concentration or response with no residual error
  vector[n_total] dv;    // concentration or response with residual error
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  
  vector[n_subjects] E0;
  vector[n_subjects] KE0;
  vector[n_subjects] EC50;
  
  {
  
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA});
    
    matrix[n_subjects, n_random] eta;   
    matrix[n_subjects, n_random] theta; 
    
    row_vector[n_random_pd] typical_values_pd = to_row_vector({TVE0, TVKE0, 
                                                               TVEC50});
    
    matrix[n_subjects, n_random_pd] eta_pd;   
    matrix[n_subjects, n_random_pd] theta_pd; 
  
    matrix[n_total, n_cmt + n_cmt_pd] x_ipred;
    
    for(i in 1:n_subjects){
      eta[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                           diag_pre_multiply(omega, L))';
                                           
      eta_pd[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random_pd),
                                              diag_pre_multiply(omega_pd, L_pd))';                                           
    }
    theta = (rep_matrix(typical_values, n_subjects) .* exp(eta));
    theta_pd = (rep_matrix(typical_values_pd, n_subjects) .* exp(eta_pd));

    CL = col(theta, 1);
    VC = col(theta, 2);
    KA = col(theta, 3);
    
    E0 = col(theta_pd, 1);
    KE0 = col(theta_pd, 2);
    EC50 = col(theta_pd, 3);
    
    for(j in 1:n_subjects){
      
      if(solver == 1){
      
        matrix[n_cmt + n_cmt_pd, n_cmt + n_cmt_pd] K = 
                              rep_matrix(0, n_cmt + n_cmt_pd, n_cmt + n_cmt_pd);
  
        K[1, 1] = -KA[j];
        K[2, 1] = KA[j];
        K[2, 2] = -CL[j]/VC[j];
        K[3, 2] = KE0[j];
        K[3, 3] = -KE0[j];
        
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
        
      }
      else if(solver == 2){
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(depot_1cmt_effcmt_ode,
                         n_cmt + n_cmt_pd,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j], KA[j], KE0[j]})';
                           
      }else{
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_onecpt_rk45(depot_1cmt_effcmt_ode_coupled,
                                n_cmt_pd,
                                time[subj_start[j]:subj_end[j]],
                                amt[subj_start[j]:subj_end[j]],
                                rate[subj_start[j]:subj_end[j]],
                                ii[subj_start[j]:subj_end[j]],
                                evid[subj_start[j]:subj_end[j]],
                                cmt[subj_start[j]:subj_end[j]],
                                addl[subj_start[j]:subj_end[j]],
                                ss[subj_start[j]:subj_end[j]],
                                {CL[j], VC[j], KA[j], KE0[j]},
                                bioav, tlag)';
                         
      }

      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){
          ipred[k] = x_ipred[k, 2] / VC[j];
        }else if(cmt[k] == 3){
          real conc_effcmt = x_ipred[k, 3]/VC[j];
          real inh = (emax*pow(conc_effcmt, hill))/(pow(EC50[j], hill) + 
                                                    pow(conc_effcmt, hill));
          ipred[k] = E0[j]*(1 - inh);
        }
      }
    
    }

    for(i in 1:n_total){
      if(ipred[i] == 0){
         dv[i] = 0;
      }else{
        
        if(cmt[i] == 2 || cmt[i] == 3){
          real ipred_tmp = ipred[i];
          real sigma_tmp = cmt[i] == 2 ? 
            sqrt(square(ipred_tmp) * Sigma[1, 1] + Sigma[2, 2] + 2*ipred_tmp*Sigma[2, 1]) :
            sqrt(square(ipred_tmp) * Sigma_pd[1, 1] + Sigma_pd[2, 2] + 2*ipred_tmp*Sigma_pd[2, 1]);
          
          dv[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
          
        }
      }
    }
  }
}



