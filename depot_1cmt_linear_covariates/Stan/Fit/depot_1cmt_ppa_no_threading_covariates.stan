// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// IIV on CL, VC, and Ka (full covariance matrix)
// proportional plus additive error - DV = IPRED*(1 + eps_p) + eps_a
// Matrix-exponential solution using Torsten (the matrix-exponential seems to be
//   faster than the analytical solution for this model)
// Deals with BLOQ values by the "CDF trick" (M4)
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   truncates the likelihood below at 0
// For PPC, it generates values from a normal that is truncated below at 0
// Covariates: 
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
  
}
data{
  
  int n_subjects;
  int n_total;
  int n_obs;
  array[n_obs] int i_obs;
  array[n_total] int ID;
  array[n_total] real amt;
  array[n_total] int cmt;
  array[n_total] int evid;
  array[n_total] real rate;
  array[n_total] real ii;
  array[n_total] int addl;
  array[n_total] int ss;
  array[n_total] real time;
  vector<lower = 0>[n_total] dv;
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;
  vector[n_total] lloq;
  array[n_total] int bloq;
  
  array[n_subjects] real<lower = 0> wt;              // baseline bodyweight (kg)
  array[n_subjects] int<lower = 0, upper = 1> cmppi; // cmppi
  array[n_subjects] real<lower = 0> egfr;            // eGFR
  
  real<lower = 0> location_tvcl;  // Prior Location parameter for CL
  real<lower = 0> location_tvvc;  // Prior Location parameter for VC
  real<lower = 0> location_tvka;  // Prior Location parameter for KA
  
  real<lower = 0> scale_tvcl;     // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;     // Prior Scale parameter for VC
  real<lower = 0> scale_tvka;     // Prior Scale parameter for KA
  
  real<lower = 0> scale_omega_cl; // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc; // Prior scale parameter for omega_vc
  real<lower = 0> scale_omega_ka; // Prior scale parameter for omega_ka
  
  real<lower = 0> lkj_df_omega;   // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma_p;  // Prior Scale parameter for proportional error
  real<lower = 0> scale_sigma_a;  // Prior Scale parameter for additive error
  
  real<lower = 0> lkj_df_sigma;   // Prior degrees of freedom for sigma cor mat
  
  int<lower = 0, upper = 1> prior_only; // Want to simulate from the prior?
 
}
transformed data{ 
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 3;                    // Number of random effects
  int n_cmt = 2;                       // Number of states in the ODEs
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc, 
                                      scale_omega_ka}; 
                                      
  array[2] real scale_sigma = {scale_sigma_p, scale_sigma_a};
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
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
transformed parameters{
  
  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_vc;
  vector[n_subjects] eta_ka;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  vector[n_subjects] KE;
  
  real<lower = 0> sigma_p = sigma[1];
  real<lower = 0> sigma_a = sigma[2];
  
  real<lower = 0> sigma_sq_p = square(sigma_p);
  real<lower = 0> sigma_sq_a = square(sigma_a);
  
  real cor_p_a;
  real sigma_p_a;
  
  vector[n_obs] ipred;

  {
  
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA});

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
    
    vector[n_total] dv_ipred;
    matrix[n_total, 2] x_ipred;
    
    eta_cl = col(eta, 1);
    eta_vc = col(eta, 2);
    eta_ka = col(eta, 3);
    
    cor_p_a = R_Sigma[1, 2];
    sigma_p_a = Sigma[1, 2];
    
    for(j in 1:n_subjects){
      
      matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
      
      real wt_over_70 = wt[j]/70;
      real wt_adjustment_cl = wt_over_70^theta_cl_wt;
      real wt_adjustment_vc = wt_over_70^theta_vc_wt;
      real cmppi_adjustment_ka = exp(theta_ka_cmppi*cmppi[j]);
      real egfr_adjustment_cl = (egfr[j]/90)^theta_cl_egfr;
      
      row_vector[n_random] theta_j = theta[j]; // access the parameters for subject j
      CL[j] = theta_j[1] * wt_adjustment_cl * egfr_adjustment_cl;
      VC[j] = theta_j[2] * wt_adjustment_vc;
      KA[j] = theta_j[3] * cmppi_adjustment_ka;
      KE[j] = CL[j]/VC[j];  
      
      K[1, 1] = -KA[j];
      K[2, 1] = KA[j];
      K[2, 2] = -CL[j]/VC[j];
      
      x_ipred[subj_start[j]:subj_end[j], ] =
        pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         K, bioav, tlag)';
                           
      dv_ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
    }
  
    ipred = dv_ipred[i_obs];
  
  }
  
}
model{ 
  
  // Priors
  TVCL ~ lognormal(log(location_tvcl), scale_tvcl);
  TVVC ~ lognormal(log(location_tvvc), scale_tvvc);
  TVKA ~ lognormal(log(location_tvka), scale_tvka) T[TVCL/TVVC, ];
  
  theta_cl_wt ~ std_normal();
  theta_vc_wt ~ std_normal();
  theta_ka_cmppi ~ std_normal();
  theta_cl_egfr ~ std_normal();

  omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma ~ normal(0, scale_sigma);
  L_Sigma ~ lkj_corr_cholesky(lkj_df_sigma);
  
  to_vector(Z) ~ std_normal();
  
  // Likelihood
  if(prior_only == 0){
    for(i in 1:n_obs){
      real ipred_tmp = ipred[i];
      real sigma_tmp = sqrt(square(ipred_tmp) * sigma_sq_p + sigma_sq_a + 
                            2*ipred_tmp*sigma_p_a);
      if(bloq_obs[i] == 1){
        target += log_diff_exp(normal_lcdf(lloq_obs[i] | ipred_tmp, sigma_tmp),
                               normal_lcdf(0.0 | ipred_tmp, sigma_tmp)) -
                   normal_lccdf(0.0 | ipred_tmp, sigma_tmp); 
      }else{
        target += normal_lpdf(dv_obs[i] | ipred_tmp, sigma_tmp) -
                  normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
      }
    }
  }
}
generated quantities{

  real<lower = 0> omega_cl = omega[1];
  real<lower = 0> omega_vc = omega[2];
  real<lower = 0> omega_ka = omega[3];

  real<lower = 0> omega_sq_cl = square(omega_cl);
  real<lower = 0> omega_sq_vc = square(omega_vc);
  real<lower = 0> omega_sq_ka = square(omega_ka);

  real cor_cl_vc;
  real cor_cl_ka;
  real cor_vc_ka;
  real omega_cl_vc;
  real omega_cl_ka;
  real omega_vc_ka;

  vector[n_obs] pred;
  vector[n_obs] dv_ppc;
  vector[n_obs] log_lik;
  vector[n_obs] res;
  vector[n_obs] wres;
  vector[n_obs] ires;
  vector[n_obs] iwres;
 
  {

    matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
    matrix[n_random, n_random] Omega = quad_form_diag(R, omega);

    vector[n_total] dv_pred;
    matrix[n_total, n_cmt] x_pred;

    cor_cl_vc = R[1, 2];
    cor_cl_ka = R[1, 3];
    cor_vc_ka = R[2, 3];

    omega_cl_vc = Omega[1, 2];
    omega_cl_ka = Omega[1, 3];
    omega_vc_ka = Omega[2, 3];

    for(j in 1:n_subjects){
      
      real wt_over_70 = wt[j]/70;
      real wt_adjustment_cl = wt_over_70^theta_cl_wt;
      real wt_adjustment_vc = wt_over_70^theta_vc_wt;
      real cmppi_adjustment_ka = exp(theta_ka_cmppi*cmppi[j]);
      real egfr_adjustment_cl = (egfr[j]/90)^theta_cl_egfr;
      
      real cl_p = TVCL * wt_adjustment_cl * egfr_adjustment_cl;
      real vc_p = TVVC * wt_adjustment_vc;
      real ka_p = TVKA * cmppi_adjustment_ka;

      matrix[n_cmt, n_cmt] K_tv = rep_matrix(0, n_cmt, n_cmt);
                           
      K_tv[1, 1] = -ka_p;
      K_tv[2, 1] = ka_p;
      K_tv[2, 2] = -cl_p/vc_p;

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         K_tv, bioav, tlag)';
      
      dv_pred[subj_start[j]:subj_end[j]] =
        x_pred[subj_start[j]:subj_end[j], 2] ./ vc_p;
      
    }

    pred = dv_pred[i_obs];

  }

  res = dv_obs - pred;
  ires = dv_obs - ipred;

  for(i in 1:n_obs){
    real ipred_tmp = ipred[i];
    real sigma_tmp = sqrt(square(ipred_tmp) * sigma_sq_p + sigma_sq_a + 
                          2*ipred_tmp*sigma_p_a);
    dv_ppc[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
    if(bloq_obs[i] == 1){
      // log_lik[i] = log(normal_cdf(lloq_obs[i] | ipred_tmp, sigma_tmp) -
      //                  normal_cdf(0.0 | ipred_tmp, sigma_tmp)) -
      //              normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
      log_lik[i] = log_diff_exp(normal_lcdf(lloq_obs[i] | ipred_tmp, sigma_tmp),
                                normal_lcdf(0.0 | ipred_tmp, sigma_tmp)) -
                   normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
    }else{
      log_lik[i] = normal_lpdf(dv_obs[i] | ipred_tmp, sigma_tmp) -
                   normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
    }
    wres[i] = res[i]/sigma_tmp;
    iwres[i] = ires[i]/sigma_tmp;
  }
  
}

