// First Order Absorption (oral/subcutaneous)
// Two-compartment PK Model
// Zero-order distributive delay (like an infusion into the gut) to bring about 
//   a delayed absorption
// IIV on CL, VC, Q, VP, KA, DUR (full covariance matrix)
// proportional error - DV = IPRED(1 + eps_p)
// Analytical solution using Torsten 
// Deals with BLOQ values by the "CDF trick" (M4)
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   truncates the likelihood below at 0
// For PPC, it generates values from a normal that is truncated below at 0

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
  array[n_total] real ii;
  array[n_total] int addl;
  array[n_total] int ss;
  array[n_total] real time;
  vector<lower = 0>[n_total] dv;
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;
  vector[n_total] lloq;
  array[n_total] int bloq;
  
  real<lower = 0> location_tvcl;   // Prior Location parameter for CL
  real<lower = 0> location_tvvc;   // Prior Location parameter for VC
  real<lower = 0> location_tvq;    // Prior Location parameter for Q
  real<lower = 0> location_tvvp;   // Prior Location parameter for VP
  real<lower = 0> location_tvka;   // Prior Location parameter for KA
  real<lower = 0> location_tvdur;  // Prior Location parameter for DUR
  
  real<lower = 0> scale_tvcl;      // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;      // Prior Scale parameter for VC
  real<lower = 0> scale_tvq;       // Prior Scale parameter for Q
  real<lower = 0> scale_tvvp;      // Prior Scale parameter for VP
  real<lower = 0> scale_tvka;      // Prior Scale parameter for KA
  real<lower = 0> scale_tvdur;     // Prior Scale parameter for DUR
  
  real<lower = 0> scale_omega_cl;  // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc;  // Prior scale parameter for omega_vc
  real<lower = 0> scale_omega_q;  // Prior scale parameter for omega_q
  real<lower = 0> scale_omega_vp;  // Prior scale parameter for omega_vp
  real<lower = 0> scale_omega_ka;  // Prior scale parameter for omega_ka
  real<lower = 0> scale_omega_dur; // Prior scale parameter for omega_dur
  
  real<lower = 0> lkj_df_omega;    // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma_p;   // Prior Scale parameter for proportional error
  
  int<lower = 0, upper = 1> prior_only; // Want to simulate from the prior?
  int<lower = 0, upper = prior_only> no_gq_predictions; // Leave out PREDS and IPREDS in 
                                                        // generated quantities. Useful
                                                        // for simulating prior parameters
                                                        // but don't want prior predictions
 
}
transformed data{ 
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 6;                    // Number of random effects
  int n_cmt = 3;                       // Number of states in the ODEs
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc,
                                      scale_omega_q, scale_omega_vp,
                                      scale_omega_ka, scale_omega_dur}; 
  
  array[n_subjects] int seq_subj = linspaced_int_array(n_subjects, 1, n_subjects); // reduce_sum over subjects
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVQ;       
  real<lower = 0> TVVP;
  real<lower = 0.5*(TVCL/TVVC + TVQ/TVVC + TVQ/TVVP +
        sqrt((TVCL/TVVC + TVQ/TVVC + TVQ/TVVP)^2 - 4*TVCL/TVVC*TVQ/TVVP))> TVKA;
  real<lower = 0> TVDUR; 
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
}
transformed parameters{
  
  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_vc;
  vector[n_subjects] eta_q;
  vector[n_subjects] eta_vp;
  vector[n_subjects] eta_ka;
  vector[n_subjects] eta_dur;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;
  vector[n_subjects] KA;
  vector[n_subjects] DUR;
  
  vector[n_obs] ipred;

  {
  
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP, 
                                                         TVKA, TVDUR});

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
    
    vector[n_total] dv_ipred;
    matrix[n_total, n_cmt] x_ipred; 
    
    eta_cl = col(eta, 1);
    eta_vc = col(eta, 2);
    eta_q = col(eta, 3);
    eta_vp = col(eta, 4);
    eta_ka = col(eta, 5);
    eta_dur = col(eta, 6);
    CL = col(theta, 1);
    VC = col(theta, 2);
    Q = col(theta, 3);
    VP = col(theta, 4);
    KA = col(theta, 5);
    DUR = col(theta, 6);
    
    for(j in 1:n_subjects){
      
      array[subj_end[j] - subj_start[j] + 1] real rate = 
              to_array_1d(to_vector(amt[subj_start[j]:subj_end[j]]) ./ DUR[j]);
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate,
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], Q[j], VC[j], VP[j], KA[j]},
                         bioav, tlag)';  
                           
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
  TVQ ~ lognormal(log(location_tvq), scale_tvq);
  TVVP ~ lognormal(log(location_tvvp), scale_tvvp);
  TVKA ~ lognormal(log(location_tvka), scale_tvka) T[0.5*(TVCL/TVVC + TVQ/TVVC + TVQ/TVVP +
          sqrt((TVCL/TVVC + TVQ/TVVC + TVQ/TVVP)^2 - 4*TVCL/TVVC*TVQ/TVVP)), ];
  TVDUR ~ lognormal(log(location_tvdur), scale_tvdur);

  omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma_p ~ normal(0, scale_sigma_p);
  
  to_vector(Z) ~ std_normal();
  
  // Likelihood
  if(prior_only == 0){
    for(i in 1:n_obs){
      real sigma_tmp = ipred[i]*sigma_p;
      if(bloq_obs[i] == 1){
        target += log_diff_exp(normal_lcdf(lloq_obs[i] | ipred[i], sigma_tmp),
                               normal_lcdf(0.0 | ipred[i], sigma_tmp)) -
                   normal_lccdf(0.0 | ipred[i], sigma_tmp); 
      }else{
        target += normal_lpdf(dv_obs[i] | ipred[i], sigma_tmp) -
                  normal_lccdf(0.0 | ipred[i], sigma_tmp);
      }
    }
  }
}
generated quantities{
  
  real<lower = 0> sigma_sq_p = square(sigma_p);

  real<lower = 0> omega_cl = omega[1];
  real<lower = 0> omega_vc = omega[2];
  real<lower = 0> omega_q = omega[3];
  real<lower = 0> omega_vp = omega[4];
  real<lower = 0> omega_ka = omega[5];
  real<lower = 0> omega_dur = omega[6];

  real<lower = 0> omega_sq_cl = square(omega_cl);
  real<lower = 0> omega_sq_vc = square(omega_vc);
  real<lower = 0> omega_sq_q = square(omega_q);
  real<lower = 0> omega_sq_vp = square(omega_vp);
  real<lower = 0> omega_sq_ka = square(omega_ka);
  real<lower = 0> omega_sq_dur = square(omega_dur);

  real cor_cl_vc;
  real cor_cl_q;
  real cor_cl_vp;
  real cor_cl_ka;
  real cor_cl_dur;
  real cor_vc_q;
  real cor_vc_vp;
  real cor_vc_ka;
  real cor_vc_dur;
  real cor_q_vp;
  real cor_q_ka;
  real cor_q_dur;
  real cor_vp_ka;
  real cor_vp_dur;
  real cor_ka_dur;
  
  real omega_cl_vc;
  real omega_cl_q;
  real omega_cl_vp;
  real omega_cl_ka;
  real omega_cl_dur;
  real omega_vc_q;
  real omega_vc_vp;
  real omega_vc_ka;
  real omega_vc_dur;
  real omega_q_vp;
  real omega_q_ka;
  real omega_q_dur;
  real omega_vp_ka;
  real omega_vp_dur;
  real omega_ka_dur;

  vector[no_gq_predictions ? 0 : n_obs] pred;
  vector[no_gq_predictions ? 0 : n_obs] epred_stan;
  vector[no_gq_predictions ? 0 : n_obs] epred;
  vector[no_gq_predictions ? 0 : n_obs] dv_ppc;
  vector[no_gq_predictions ? 0 : n_obs] log_lik;
  vector[no_gq_predictions ? 0 : n_obs] iwres;
 
  {
    
    matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
    matrix[n_random, n_random] Omega = quad_form_diag(R, omega);

    cor_cl_vc = R[1, 2];
    cor_cl_q = R[1, 3];
    cor_cl_vp = R[1, 4];
    cor_cl_ka = R[1, 5];
    cor_cl_dur = R[1, 6];
    cor_vc_q = R[2, 3];
    cor_vc_vp = R[2, 4];
    cor_vc_ka = R[2, 5];
    cor_vc_dur = R[2, 6];
    cor_q_vp = R[3, 4];
    cor_q_ka = R[3, 5];
    cor_q_dur = R[3, 6];
    cor_vp_ka = R[4, 5];
    cor_vp_dur = R[4, 6];
    cor_ka_dur = R[5, 6];

    omega_cl_vc = Omega[1, 2];
    omega_cl_q = Omega[1, 3];
    omega_cl_vp = Omega[1, 4];
    omega_cl_ka = Omega[1, 5];
    omega_cl_dur = Omega[1, 6];
    omega_vc_q = Omega[2, 3];
    omega_vc_vp = Omega[2, 4];
    omega_vc_ka = Omega[2, 5];
    omega_vc_dur = Omega[2, 6];
    omega_q_vp = Omega[3, 4];
    omega_q_ka = Omega[3, 5];
    omega_q_dur = Omega[3, 6];
    omega_vp_ka = Omega[4, 5];
    omega_vp_dur = Omega[4, 6];
    omega_ka_dur = Omega[5, 6];
    
  }
  
  if(no_gq_predictions == 0){
    
    vector[n_subjects] CL_new;
    vector[n_subjects] VC_new;
    vector[n_subjects] Q_new;
    vector[n_subjects] VP_new;
    vector[n_subjects] KA_new;
    vector[n_subjects] DUR_new;
    
    vector[n_total] dv_pred;
    matrix[n_total, n_cmt] x_pred;
    vector[n_total] dv_epred;
    matrix[n_total, n_cmt] x_epred;
    
    matrix[n_subjects, n_random] eta_new;
    matrix[n_subjects, n_random] theta_new;
    
    array[n_total] real rate_new;
    array[n_total] real rate_p;
    
    for(i in 1:n_subjects){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = 
      (rep_matrix(to_row_vector({TVCL, TVVC, TVQ, TVVP, TVKA, TVDUR}), n_subjects) .* 
                                                                  exp(eta_new));

    for(j in 1:n_subjects){
    
      row_vector[n_random] theta_j_new = theta_new[j]; // access the parameters for subject j's epred
      
      real cl_p = TVCL;
      real vc_p = TVVC;
      real q_p = TVQ;
      real vp_p = TVVP;
      real ka_p = TVKA;
      real dur_p = TVDUR;
      
      CL_new[j] = theta_j_new[1];
      VC_new[j] = theta_j_new[2];
      Q_new[j] = theta_j_new[3];
      VP_new[j] = theta_j_new[4];
      KA_new[j] = theta_j_new[5];
      DUR_new[j] = theta_j_new[6];
              
      rate_new[subj_start[j]:subj_end[j]] =  
           to_array_1d(to_vector(amt[subj_start[j]:subj_end[j]]) ./ DUR_new[j]);
           
      rate_p[subj_start[j]:subj_end[j]] =  
           to_array_1d(to_vector(amt[subj_start[j]:subj_end[j]]) ./ dur_p);
      
      x_epred[subj_start[j]:subj_end[j], ] =
        pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate_new[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL_new[j], Q_new[j], VC_new[j], VP_new[j], KA_new[j]},
                         bioav, tlag)';
                         
      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate_p[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {cl_p, q_p, vc_p, vp_p, ka_p},
                         bioav, tlag)';
        
      dv_epred[subj_start[j]:subj_end[j]] =
        x_epred[subj_start[j]:subj_end[j], 2] ./ VC_new[j];
      
      dv_pred[subj_start[j]:subj_end[j]] =
        x_pred[subj_start[j]:subj_end[j], 2] ./ vc_p;
      
    }

    pred = dv_pred[i_obs];
    epred_stan = dv_epred[i_obs];

    for(i in 1:n_obs){
      real ipred_tmp = ipred[i];
      real sigma_tmp = ipred_tmp*sigma_p;
      real epred_tmp = epred_stan[i];
      real sigma_tmp_e = epred_tmp*sigma_p;
      dv_ppc[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
      epred[i] = normal_lb_rng(epred_tmp, sigma_tmp_e, 0.0);
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
      iwres[i] = (dv_obs[i] - ipred_tmp)/sigma_tmp;
    }
  }
}
