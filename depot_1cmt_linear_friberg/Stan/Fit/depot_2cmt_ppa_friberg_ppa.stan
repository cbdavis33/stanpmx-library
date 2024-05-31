// First Order Absorption (oral/subcutaneous)
// Two-compartment PK Model, Friberg-Karlsson neutropenia PD model
// IIV on CL, VC, Q, VP, and KA (full covariance matrix)
// IIV on MTT, CIRC0, GAMMA, ALPHA (full covariance matrix) for PD
// proportional plus additive error on PK - DV = IPRED*(1 + eps_p) + eps_a
// proportional plus additive error on PD - DV = IPRED*(1 + eps_p_pd) + eps_a_pd
// General ODE solution using Torsten (not sure why, but it's faster than 
//   coupled). rk45, since it's faster than bdf
// Implements threading for within-chain parallelization 
// Deals with BLOQ values by the "CDF trick" (M4)
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   truncates the likelihood below at 0
// For PPC, it generates values from a normal that is truncated below at 0

functions{

  array[] int sequence(int start, int end) { 
    array[end - start + 1] int seq;
    for (n in 1:num_elements(seq)) {
      seq[n] = n + start - 1;
    }
    return seq; 
  } 
  
  int num_between(int lb, int ub, array[] int y){
    
    int n = 0;
    for(i in 1:num_elements(y)){
      if(y[i] >= lb && y[i] <= ub)
         n = n + 1;
    }
    return n;
    
  }
  
  array[] int find_between(int lb, int ub, array[] int y) {
    // vector[num_between(lb, ub, y)] result;
    array[num_between(lb, ub, y)] int result;
    int n = 1;
    for (i in 1:num_elements(y)) {
      if (y[i] >= lb && y[i] <= ub) {
        result[n] = y[i];
        n = n + 1;
      }
    }
    return result;
  }
  
  vector find_between_vec(int lb, int ub, array[] int idx, vector y) {
    
    vector[num_between(lb, ub, idx)] result;
    int n = 1;
    if(num_elements(idx) != num_elements(y)) reject("illegal input");
    for (i in 1:rows(y)) {
      if (idx[i] >= lb && idx[i] <= ub) {
        result[n] = y[i];
        n = n + 1;
      }
    }
    return result;
  }
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y; 

  }
  
  vector depot_2cmt_friberg_ode(real t, vector y, array[] real params, 
                                array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real q = params[3];
    real vp = params[4];
    real ka = params[5];
    real mtt = params[6];
    real circ_0 = params[7];
    real gamma = params[8];
    real alpha = params[9];
    
    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;
    real k_tr = 4/mtt; // k_tr = (n_tr + 1)/mtt    
    
    real conc = y[2]/vc;
    
    real e_drug = fmin(1.0, alpha*conc); // Maybe reparameterize this so no more fmin?
    real prol = y[4] + circ_0;
    real transit_1 = y[5] + circ_0; 
    real transit_2 = y[6] + circ_0;
    real transit_3 = y[7] + circ_0;
    real circ = y[8] + circ_0; // fmax(machine_precision(), y[8] + circ_0)
    
    vector[8] dydt;
    
    dydt[1] = -ka*y[1];                               // depot
    dydt[2] = ka*y[1] - (ke + k_cp)*y[2] + k_pc*y[3]; // central
    dydt[3] = k_cp*y[2] - k_pc*y[3];                  // peripheral
    dydt[4] = k_tr*prol*((1 - e_drug)*(circ_0/circ)^gamma - 1);  // proliferative cells
    dydt[5] = k_tr*(prol - transit_1);                // transit 1
    dydt[6] = k_tr*(transit_1 - transit_2);           // transit 2
    dydt[7] = k_tr*(transit_2 - transit_3);           // transit 3
    dydt[8] = k_tr*(transit_3 - circ);                // circulating blood cells
    
    return dydt;
  }
  
  vector depot_2cmt_friberg_ode_coupled(real t, vector y, vector y_pk,
                                        array[] real params, array[] real x_r, 
                                        array[] int x_i){
    
    real vc = params[3];
    real mtt = params[6];
    real circ_0 = params[7];
    real gamma = params[8];
    real alpha = params[9];
    
    real conc = y_pk[2]/vc;
    
    real k_tr = 4/mtt; // k_tr = (n_tr + 1)/mtt    
    real e_drug = fmin(1.0, alpha*conc); // Maybe reparameterize this so no more fmin?
    real prol = y[1] + circ_0;
    real transit_1 = y[2] + circ_0; 
    real transit_2 = y[3] + circ_0;
    real transit_3 = y[4] + circ_0;
    real circ = y[5] + circ_0; 
    
    vector[5] dydt;
    
    dydt[1] = k_tr*prol*((1 - e_drug)*(circ_0/circ)^gamma - 1);  // proliferative cells
    dydt[2] = k_tr*(prol - transit_1);                // transit 1
    dydt[3] = k_tr*(transit_1 - transit_2);           // transit 2
    dydt[4] = k_tr*(transit_2 - transit_3);           // transit 3
    dydt[5] = k_tr*(transit_3 - circ);                // circulating blood cells
    
    return dydt;
  }
  
  real partial_sum_lpmf(array[] int seq_subj, int start, int end,
                        vector dv_obs, array[] int dv_obs_id, array[] int i_obs,
                        array[] real amt, array[] int cmt, array[] int evid, 
                        array[] real time, array[] real rate, array[] real ii, 
                        array[] int addl, array[] int ss,
                        array[] int subj_start, array[] int subj_end, 
                        vector CL, vector VC, vector Q, vector VP, vector KA, 
                        vector MTT, vector CIRC0, vector GAMMA, vector ALPHA,
                        real sigma_sq_p, real sigma_sq_a, real sigma_p_a,
                        real sigma_sq_p_pd, real sigma_sq_a_pd, real sigma_p_a_pd,
                        vector lloq, array[] int bloq,
                        int n_random, int n_random_pd, 
                        int n_subjects, int n_total,
                        array[] real bioav, array[] real tlag, 
                        int n_cmt, int n_cmt_pd){
                           
    real ptarget = 0;
                              
    int N = end - start + 1;    // number of subjects in this slice  
    vector[n_total] dv_ipred;   
    matrix[n_total, n_cmt + n_cmt_pd] x_ipred;
  
    int n_obs_slice = num_between(subj_start[start], subj_end[end], i_obs);
    array[n_obs_slice] int i_obs_slice = find_between(subj_start[start], 
                                                      subj_end[end], i_obs);
                                                
    vector[n_obs_slice] dv_obs_slice = find_between_vec(start, end, 
                                                        dv_obs_id, dv_obs);
    
    vector[n_obs_slice] ipred_slice;
    
    vector[n_obs_slice] lloq_slice = lloq[i_obs_slice];
    array[n_obs_slice] int bloq_slice = bloq[i_obs_slice];
    
    array[n_obs_slice] int cmt_slice = cmt[i_obs_slice];
    
    for(n in 1:N){            // loop over subjects in this slice
    
      int j = n + start - 1; // j is the ID of the current subject
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_2cmt_friberg_ode,
                       n_cmt + n_cmt_pd,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {CL[j], VC[j], Q[j], VP[j], KA[j], 
                        MTT[j], CIRC0[j], GAMMA[j], ALPHA[j]})';
          
      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){
          dv_ipred[k] = x_ipred[k, 2] / VC[j];
        }else if(cmt[k] == 4){
          dv_ipred[k] = x_ipred[k, 8] + CIRC0[j];
        }
      }
    
    }
  
    ipred_slice = dv_ipred[i_obs_slice];
    
    for(i in 1:n_obs_slice){
    
      if(cmt_slice[i] == 2 || cmt_slice[i] == 4){
        real ipred_tmp = ipred_slice[i];
        real sigma_tmp = cmt_slice[i] == 2 ? 
          sqrt(square(ipred_tmp) * sigma_sq_p + sigma_sq_a + 2*ipred_tmp*sigma_p_a) : 
          sqrt(square(ipred_tmp) * sigma_sq_p_pd + sigma_sq_a_pd + 2*ipred_tmp*sigma_p_a_pd);
    
        if(bloq_slice[i] == 1){
          ptarget += log_diff_exp(normal_lcdf(lloq_slice[i] | ipred_tmp, sigma_tmp),
                                  normal_lcdf(0.0 | ipred_tmp, sigma_tmp)) -
                     normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
        }else{
          ptarget += normal_lpdf(dv_obs_slice[i] | ipred_tmp, sigma_tmp) -
                     normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
        }
      }  
    }                                         
                              
    return ptarget;
                           
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
  
  real<lower = 0> location_tvcl;  // Prior Location parameter for CL
  real<lower = 0> location_tvvc;  // Prior Location parameter for VC
  real<lower = 0> location_tvq;   // Prior Location parameter for Q
  real<lower = 0> location_tvvp;  // Prior Location parameter for VP
  real<lower = 0> location_tvka;  // Prior Location parameter for KA
  
  real<lower = 0> scale_tvcl;     // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;     // Prior Scale parameter for VC
  real<lower = 0> scale_tvq;      // Prior Scale parameter for Q
  real<lower = 0> scale_tvvp;     // Prior Scale parameter for VP
  real<lower = 0> scale_tvka;     // Prior Scale parameter for KA
  
  real<lower = 0> scale_omega_cl; // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc; // Prior scale parameter for omega_vc
  real<lower = 0> scale_omega_q;  // Prior scale parameter for omega_q
  real<lower = 0> scale_omega_vp; // Prior scale parameter for omega_vp
  real<lower = 0> scale_omega_ka; // Prior scale parameter for omega_ka
  
  real<lower = 0> lkj_df_omega;   // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma_p;  // Prior Scale parameter for proportional error
  real<lower = 0> scale_sigma_a;  // Prior Scale parameter for additive error
  
  real<lower = 0> lkj_df_sigma;   // Prior degrees of freedom for sigma cor mat
  
  real<lower = 0> location_tvmtt;   // Prior Location parameter for MTT
  real<lower = 0> location_tvcirc0; // Prior Location parameter for CIRC0
  real<lower = 0> location_tvgamma; // Prior Location parameter for GAMMA
  real<lower = 0> location_tvalpha; // Prior Location parameter for ALPHA
  
  real<lower = 0> scale_tvmtt;      // Prior Scale parameter for MTT
  real<lower = 0> scale_tvcirc0;    // Prior Scale parameter for CIRC0
  real<lower = 0> scale_tvgamma;    // Prior Scale parameter for GAMMA
  real<lower = 0> scale_tvalpha;    // Prior Scale parameter for ALPHA
  
  real<lower = 0> scale_omega_mtt;   // Prior scale parameter for omega_mtt
  real<lower = 0> scale_omega_circ0; // Prior scale parameter for omega_circ0
  real<lower = 0> scale_omega_gamma; // Prior scale parameter for omega_gamma
  real<lower = 0> scale_omega_alpha; // Prior scale parameter for omega_alpha
  
  real<lower = 0> lkj_df_omega_pd;  // Prior degrees of freedom for omega_pd cor mat
  
  real<lower = 0> scale_sigma_p_pd;  // Prior Scale parameter for proportional error
  real<lower = 0> scale_sigma_a_pd;  // Prior Scale parameter for additive error
  
  real<lower = 0> lkj_df_sigma_pd;   // Prior degrees of freedom for sigma_pd cor mat
  
  int<lower = 0, upper = 1> prior_only; // Want to simulate from the prior?
  int<lower = 0, upper = prior_only> no_gq_predictions; // Leave out PREDS and IPREDS in 
                                                        // generated quantities. Useful
                                                        // for simulating prior parameters
                                                        // but don't want prior predictions
 
}
transformed data{ 
  
  int grainsize = 1;
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 5;
  int n_random_pd = 4;
  
  int n_cmt = 3;     // number of ODEs in PK model (depot, central, peripheral)
  int n_cmt_pd = 5;  // number of ODEs in PD system
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc, 
                                      scale_omega_q, scale_omega_vp,
                                      scale_omega_ka};
                                      
  array[2] real scale_sigma = {scale_sigma_p, scale_sigma_a};
                                      
  array[n_random_pd] real scale_omega_pd = {scale_omega_mtt, scale_omega_circ0, 
                                            scale_omega_gamma, scale_omega_alpha};
                                            
  array[2] real scale_sigma_pd = {scale_sigma_p_pd, scale_sigma_a_pd};
  
  array[n_subjects] int seq_subj = sequence(1, n_subjects); // reduce_sum over subjects
  
  array[n_cmt + n_cmt_pd] real bioav = rep_array(1.0, n_cmt + n_cmt_pd); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt + n_cmt_pd] real tlag = rep_array(0.0, n_cmt + n_cmt_pd);
  
}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVQ;       
  real<lower = 0> TVVP;
  real<lower = 0.5*(TVCL/TVVC + TVQ/TVVC + TVQ/TVVP +
    sqrt((TVCL/TVVC + TVQ/TVVC + TVQ/TVVP)^2 - 4*TVCL/TVVC*TVQ/TVVP))> TVKA;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  vector<lower = 0>[2] sigma;
  cholesky_factor_corr[2] L_Sigma;
  
  matrix[n_random, n_subjects] Z;
  
  real<lower = 0> TVMTT;       
  real<lower = 0> TVCIRC0; 
  real<lower = 0> TVGAMMA;
  real<lower = 0> TVALPHA;
  
  vector<lower = 0>[n_random_pd] omega_pd;
  cholesky_factor_corr[n_random_pd] L_pd;
  
  vector<lower = 0>[2] sigma_pd;
  cholesky_factor_corr[2] L_Sigma_pd;
  
  matrix[n_random_pd, n_subjects] Z_pd;
  
}
transformed parameters{
  
  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_vc;
  vector[n_subjects] eta_q;
  vector[n_subjects] eta_vp;
  vector[n_subjects] eta_ka;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;
  vector[n_subjects] KA;
  
  real<lower = 0> sigma_p = sigma[1];
  real<lower = 0> sigma_a = sigma[2];
  
  real<lower = 0> sigma_sq_p = square(sigma_p);
  real<lower = 0> sigma_sq_a = square(sigma_a);
  
  real cor_p_a;
  real sigma_p_a;
  
  vector[n_subjects] eta_mtt;
  vector[n_subjects] eta_circ0;
  vector[n_subjects] eta_gamma;
  vector[n_subjects] eta_alpha;
  vector[n_subjects] MTT;
  vector[n_subjects] CIRC0;
  vector[n_subjects] GAMMA;
  vector[n_subjects] ALPHA;
  
  real<lower = 0> sigma_p_pd = sigma_pd[1];
  real<lower = 0> sigma_a_pd = sigma_pd[2];
  
  real<lower = 0> sigma_sq_p_pd = square(sigma_p_pd);
  real<lower = 0> sigma_sq_a_pd = square(sigma_a_pd);
  
  real cor_p_a_pd;
  real sigma_p_a_pd;

  {
    
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP, 
                                                         TVKA});
    
    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';  
    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
    
    row_vector[n_random_pd] typical_values_pd = to_row_vector({TVMTT, TVCIRC0, 
                                                               TVGAMMA, TVALPHA});
    
    matrix[n_subjects, n_random_pd] eta_pd = 
                                      diag_pre_multiply(omega_pd, L_pd * Z_pd)';  
    matrix[n_subjects, n_random_pd] theta_pd =
                     (rep_matrix(typical_values_pd, n_subjects) .* exp(eta_pd));
                     
    matrix[2, 2] R_Sigma_pd = multiply_lower_tri_self_transpose(L_Sigma_pd);
    matrix[2, 2] Sigma_pd = quad_form_diag(R_Sigma_pd, sigma_pd);
    
    eta_cl = col(eta, 1);
    eta_vc = col(eta, 2);
    eta_q = col(eta, 3);
    eta_vp = col(eta, 4);
    eta_ka = col(eta, 5);
    CL = col(theta, 1);
    VC = col(theta, 2);
    Q = col(theta, 3);
    VP = col(theta, 4);
    KA = col(theta, 5);
  
    eta_mtt = col(eta_pd, 1);
    eta_circ0 = col(eta_pd, 2);
    eta_gamma = col(eta_pd, 3);
    eta_alpha = col(eta_pd, 4);
    MTT = col(theta_pd, 1);
    CIRC0 = col(theta_pd, 2);
    GAMMA = col(theta_pd, 3);
    ALPHA = col(theta_pd, 4);
    
    cor_p_a = R_Sigma[1, 2];
    sigma_p_a = Sigma[1, 2];
    
    cor_p_a_pd = R_Sigma_pd[1, 2];
    sigma_p_a_pd = Sigma_pd[1, 2];
    
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

  omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma ~ normal(0, scale_sigma);
  L_Sigma ~ lkj_corr_cholesky(lkj_df_sigma);
  
  to_vector(Z) ~ std_normal();
  
  TVMTT ~ lognormal(log(location_tvmtt), scale_tvmtt);
  TVCIRC0 ~ lognormal(log(location_tvcirc0), scale_tvcirc0);
  TVGAMMA ~ lognormal(log(location_tvgamma), scale_tvgamma);
  TVALPHA ~ lognormal(log(location_tvalpha), scale_tvalpha);

  omega_pd ~ normal(0, scale_omega_pd);
  L_pd ~ lkj_corr_cholesky(lkj_df_omega_pd);
  
  sigma_pd ~ normal(0, scale_sigma_pd);
  L_Sigma_pd ~ lkj_corr_cholesky(lkj_df_sigma_pd);
  
  to_vector(Z_pd) ~ std_normal();
  
  // Likelihood
  if(prior_only == 0){
    target += reduce_sum(partial_sum_lupmf, seq_subj, grainsize,
                         dv_obs, dv_obs_id, i_obs,
                         amt, cmt, evid, time, 
                         rate, ii, addl, ss, subj_start, subj_end, 
                         CL, VC, Q, VP, KA,
                         MTT, CIRC0, GAMMA, ALPHA,
                         sigma_sq_p, sigma_sq_a, sigma_p_a, 
                         sigma_sq_p_pd, sigma_sq_a_pd, sigma_p_a_pd,
                         lloq, bloq,
                         n_random, n_random_pd, n_subjects, n_total,
                         bioav, tlag, n_cmt, n_cmt_pd);
                         
  }
}
generated quantities{

  real<lower = 0> omega_cl = omega[1];
  real<lower = 0> omega_vc = omega[2];
  real<lower = 0> omega_q = omega[3];
  real<lower = 0> omega_vp = omega[4];
  real<lower = 0> omega_ka = omega[5];

  real<lower = 0> omega_sq_cl = square(omega_cl);
  real<lower = 0> omega_sq_vc = square(omega_vc);
  real<lower = 0> omega_sq_q = square(omega_q);
  real<lower = 0> omega_sq_vp = square(omega_vp);
  real<lower = 0> omega_sq_ka = square(omega_ka);

  real cor_cl_vc;
  real cor_cl_q;
  real cor_cl_vp;
  real cor_cl_ka;
  real cor_vc_q;
  real cor_vc_vp;
  real cor_vc_ka;
  real cor_q_vp;
  real cor_q_ka;
  real cor_vp_ka;
  real omega_cl_vc;
  real omega_cl_q;
  real omega_cl_vp;
  real omega_cl_ka;
  real omega_vc_q;
  real omega_vc_vp;
  real omega_vc_ka;
  real omega_q_vp;
  real omega_q_ka;
  real omega_vp_ka;
  
  real<lower = 0> omega_mtt = omega_pd[1];
  real<lower = 0> omega_circ0 = omega_pd[2];
  real<lower = 0> omega_gamma = omega_pd[3];
  real<lower = 0> omega_alpha = omega_pd[4];

  real<lower = 0> omega_sq_mtt = square(omega_mtt);
  real<lower = 0> omega_sq_circ0 = square(omega_circ0);
  real<lower = 0> omega_sq_gamma = square(omega_gamma);
  real<lower = 0> omega_sq_alpha = square(omega_alpha);

  real cor_mtt_circ0;
  real cor_mtt_gamma;
  real cor_mtt_alpha;
  real cor_circ0_gamma;
  real cor_circ0_alpha;
  real cor_gamma_alpha;
  
  real omega_mtt_circ0;
  real omega_mtt_gamma;
  real omega_mtt_alpha;
  real omega_circ0_gamma;
  real omega_circ0_alpha;
  real omega_gamma_alpha;
  
  vector[n_obs] ipred;
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
    
    matrix[n_random_pd, n_random_pd] R_pd = 
                                       multiply_lower_tri_self_transpose(L_pd);
    matrix[n_random_pd, n_random_pd] Omega_pd = quad_form_diag(R_pd, omega_pd);
    
    cor_cl_vc = R[1, 2];
    cor_cl_q = R[1, 3];
    cor_cl_vp = R[1, 4];
    cor_cl_ka = R[1, 5];
    cor_vc_q = R[2, 3];
    cor_vc_vp = R[2, 4];
    cor_vc_ka = R[2, 5];
    cor_q_vp = R[3, 4];
    cor_q_ka = R[3, 5];
    cor_vp_ka = R[4, 5];

    omega_cl_vc = Omega[1, 2];
    omega_cl_q = Omega[1, 3];
    omega_cl_vp = Omega[1, 4];
    omega_cl_ka = Omega[1, 5];
    omega_vc_q = Omega[2, 3];
    omega_vc_vp = Omega[2, 4];
    omega_vc_ka = Omega[2, 5];
    omega_q_vp = Omega[3, 4];
    omega_q_ka = Omega[3, 5];
    omega_vp_ka = Omega[4, 5];
    
    cor_mtt_circ0 = R_pd[1, 2];
    cor_mtt_gamma = R_pd[1, 3];
    cor_mtt_alpha = R_pd[1, 4];
    cor_circ0_gamma = R_pd[2, 3];
    cor_circ0_alpha = R_pd[2, 4];
    cor_gamma_alpha = R_pd[3, 4];

    omega_mtt_circ0 = Omega_pd[1, 2];
    omega_mtt_gamma = Omega_pd[1, 3];
    omega_mtt_alpha = Omega_pd[1, 4];
    omega_circ0_gamma = Omega_pd[2, 3];
    omega_circ0_alpha = Omega_pd[2, 4];
    omega_gamma_alpha = Omega_pd[3, 4];
    
  }

  if(no_gq_predictions == 0){
    
    vector[n_total] dv_pred;
    matrix[n_total, n_cmt + n_cmt_pd] x_pred;
    vector[n_total] dv_ipred;
    matrix[n_total, n_cmt + n_cmt_pd] x_ipred;
    
    for(j in 1:n_subjects){
        
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_2cmt_friberg_ode,
                       n_cmt + n_cmt_pd,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {CL[j], VC[j], Q[j], VP[j], KA[j], 
                        MTT[j], CIRC0[j], GAMMA[j], ALPHA[j]})';
                        
      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_2cmt_friberg_ode,
                       n_cmt + n_cmt_pd,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {TVCL, TVVC, TVQ, TVVP, TVKA, 
                        TVMTT, TVCIRC0, TVGAMMA, TVALPHA})';
          
      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){
          dv_ipred[k] = x_ipred[k, 2] / VC[j];
          dv_pred[k] = x_pred[k, 2] / TVVC;
        }else if(cmt[k] == 4){
          dv_ipred[k] = x_ipred[k, 8] + CIRC0[j];
          dv_pred[k] = x_pred[k, 8] + TVCIRC0;
        }
      }
    }

    pred = dv_pred[i_obs];
    ipred = dv_ipred[i_obs];
    
    res = dv_obs - pred;
    ires = dv_obs - ipred;
    
    for(i in 1:n_obs){
    
      if(cmt[i_obs[i]] == 2 || cmt[i_obs[i]] == 4){
        real ipred_tmp = ipred[i];
        real sigma_tmp = cmt[i_obs[i]] == 2 ? 
          sqrt(square(ipred_tmp) * sigma_sq_p + sigma_sq_a + 2*ipred_tmp*sigma_p_a) : 
          sqrt(square(ipred_tmp) * sigma_sq_p_pd + sigma_sq_a_pd + 
                            2*ipred_tmp*sigma_p_a_pd);
    
        dv_ppc[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
    
        if(bloq_obs[i] == 1){
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
  }
}

