// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model with Indirect Response 2 PD Model
// Inhibitory Effect Compartment PD Model
// IIV on CL, VC, and KA (full covariance matrix)
// IIV on E0, KE0, and EC50 (full covariance matrix) for PD. EMAX and HILL are 
//   fixed to be 1
// For PD - IPRED = E0*(1 - emax*conc_effcmt^hill/(ec50^hill + conc_effcmt^hill))
// exponential error on PK - DV = IPRED*exp(eps)
// exponential error on PD - DV = IPRED*exp(eps_pd)
// Matrix-exponential solution using Torsten
// Implements threading for within-chain parallelization 
// Deals with BLOQ values by the M3 method (M3 and M4 are equivalent with this
//   error model)

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
  
  real partial_sum_lpmf(array[] int seq_subj, int start, int end,
                        vector dv_obs, array[] int dv_obs_id, array[] int i_obs,
                        array[] real amt, array[] int cmt, array[] int evid, 
                        array[] real time, array[] real rate, array[] real ii, 
                        array[] int addl, array[] int ss,
                        array[] int subj_start, array[] int subj_end, 
                        vector CL, vector VC, vector KA, 
                        vector E0, vector KE0, vector EC50, 
                        vector EMAX, vector HILL,
                        real sigma, real sigma_pd,
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
                      
      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){
          dv_ipred[k] = x_ipred[k, 2] / VC[j];
        }else if(cmt[k] == 3){
          real conc_effcmt = x_ipred[k, 3]/VC[j];
          real inh = (EMAX[j]*pow(conc_effcmt, HILL[j]))/(pow(EC50[j], HILL[j]) + 
                                                    pow(conc_effcmt, HILL[j]));
          dv_ipred[k] = E0[j]*(1 - inh);
        }
      }
    
    }
  
    ipred_slice = dv_ipred[i_obs_slice];
    
    for(i in 1:n_obs_slice){
      
      real log_ipred_tmp = log(ipred_slice[i]);
      real sigma_tmp = cmt_slice[i] == 2 ? sigma : sigma_pd;
      
      if(cmt_slice[i] == 2 || cmt_slice[i] == 3){
        
        if(bloq_slice[i] == 1){
          ptarget += lognormal_lcdf(lloq_slice[i] | log_ipred_tmp, sigma_tmp); 
        }else{
          ptarget += lognormal_lpdf(dv_obs_slice[i] | log_ipred_tmp, sigma_tmp);
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
  real<lower = 0> location_tvka;  // Prior Location parameter for KA
  
  real<lower = 0> scale_tvcl;     // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;     // Prior Scale parameter for VC
  real<lower = 0> scale_tvka;     // Prior Scale parameter for KA
  
  real<lower = 0> scale_omega_cl; // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc; // Prior scale parameter for omega_vc
  real<lower = 0> scale_omega_ka; // Prior scale parameter for omega_ka
  
  real<lower = 0> lkj_df_omega;   // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma;    // Prior Scale parameter for exponential error
  
  real<lower = 0> location_tve0;   // Prior Location parameter for E0
  real<lower = 0> location_tvke0;  // Prior Location parameter for KE0
  real<lower = 0> location_tvec50; // Prior Location parameter for EC50
  
  real<lower = 0> scale_tve0;      // Prior Scale parameter for E0
  real<lower = 0> scale_tvke0;     // Prior Scale parameter for KE0
  real<lower = 0> scale_tvec50;    // Prior Scale parameter for EC50
  
  real<lower = 0> scale_omega_e0;   // Prior scale parameter for omega_e0
  real<lower = 0> scale_omega_ke0;  // Prior scale parameter for omega_ke0
  real<lower = 0> scale_omega_ec50; // Prior scale parameter for omega_ec50
  
  real<lower = 0> lkj_df_omega_pd;  // Prior degrees of freedom for omega_pd cor mat
  
  real<lower = 0> scale_sigma_pd;  // Prior Scale parameter for exponential error for PD
  
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
  
  int n_random = 3;
  int n_random_pd = 3;
  
  int n_cmt = 2;     // number of ODEs in PK model (depot, central, peripheral)
  int n_cmt_pd = 1;  // number of ODEs in PD system
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc, 
                                      scale_omega_ka};
                                      
  array[n_random_pd] real scale_omega_pd = {scale_omega_e0, scale_omega_ke0, 
                                            scale_omega_ec50};
  
  array[n_subjects] int seq_subj = sequence(1, n_subjects); // reduce_sum over subjects
  
  array[n_cmt + n_cmt_pd] real bioav = rep_array(1.0, n_cmt + n_cmt_pd); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt + n_cmt_pd] real tlag = rep_array(0.0, n_cmt + n_cmt_pd);
  
  real TVEMAX = 1.0;
  real TVHILL = 1.0;
  vector<lower = 0>[n_subjects] EMAX = rep_vector(TVEMAX, n_subjects); // EMAX and HILL are both fixed to 1.0 in this model,
  vector<lower = 0>[n_subjects] HILL = rep_vector(TVHILL, n_subjects); // but it could be data or a parameter in another model.
                                                                       // Putting this here will require the least amount of 
                                                                       // changesin the code if that change is made
  
}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = TVCL/TVVC> TVKA;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma;
  
  matrix[n_random, n_subjects] Z;
  
  real<lower = 0> TVE0;       
  real<lower = 0> TVKE0; 
  real<lower = 0> TVEC50;
  
  vector<lower = 0>[n_random_pd] omega_pd;
  cholesky_factor_corr[n_random_pd] L_pd;
  
  real<lower = 0> sigma_pd;
  
  matrix[n_random_pd, n_subjects] Z_pd;
  
}
transformed parameters{
  
  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_vc;
  vector[n_subjects] eta_ka;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  
  vector[n_subjects] eta_e0;
  vector[n_subjects] eta_ke0;
  vector[n_subjects] eta_ec50;
  vector[n_subjects] E0;
  vector[n_subjects] KE0;
  vector[n_subjects] EC50;

  {
    
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA});
    
    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';  
    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
    
    row_vector[n_random_pd] typical_values_pd = to_row_vector({TVE0, TVKE0, 
                                                               TVEC50});
    
    matrix[n_subjects, n_random_pd] eta_pd = 
                                      diag_pre_multiply(omega_pd, L_pd * Z_pd)';  
    matrix[n_subjects, n_random_pd] theta_pd =
                     (rep_matrix(typical_values_pd, n_subjects) .* exp(eta_pd));
    
    eta_cl = col(eta, 1);
    eta_vc = col(eta, 2);
    eta_ka = col(eta, 3);
    CL = col(theta, 1);
    VC = col(theta, 2);
    KA = col(theta, 3);
  
    eta_e0 = col(eta_pd, 1);
    eta_ke0 = col(eta_pd, 2);
    eta_ec50 = col(eta_pd, 3);
    E0 = col(theta_pd, 1);
    KE0 = col(theta_pd, 2);
    EC50 = col(theta_pd, 3);
    
  }
  
}
model{ 
  
  // Priors
  TVCL ~ lognormal(log(location_tvcl), scale_tvcl);
  TVVC ~ lognormal(log(location_tvvc), scale_tvvc);
  TVKA ~ lognormal(log(location_tvka), scale_tvka) T[TVCL/TVVC, ];

  omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma ~ normal(0, scale_sigma);
  
  to_vector(Z) ~ std_normal();
  
  TVE0 ~ lognormal(log(location_tve0), scale_tve0);
  TVKE0 ~ lognormal(log(location_tvke0), scale_tvke0);
  TVEC50 ~ lognormal(log(location_tvec50), scale_tvec50);

  omega_pd ~ normal(0, scale_omega_pd);
  L_pd ~ lkj_corr_cholesky(lkj_df_omega_pd);
  
  sigma_pd ~ normal(0, scale_sigma_pd);
  
  to_vector(Z_pd) ~ std_normal();
  
  // Likelihood
  if(prior_only == 0){
    target += reduce_sum(partial_sum_lupmf, seq_subj, grainsize,
                         dv_obs, dv_obs_id, i_obs,
                         amt, cmt, evid, time, 
                         rate, ii, addl, ss, subj_start, subj_end, 
                         CL, VC, KA,
                         E0, KE0, EC50, EMAX, HILL,
                         sigma, sigma_pd,
                         lloq, bloq,
                         n_random, n_random_pd, n_subjects, n_total,
                         bioav, tlag, n_cmt, n_cmt_pd);
                         
  }
}
generated quantities{
  
  real<lower = 0> sigma_sq = square(sigma);
  real<lower = 0> sigma_sq_pd = square(sigma_pd);

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
  
  real<lower = 0> omega_e0 = omega_pd[1];
  real<lower = 0> omega_ke0 = omega_pd[2];
  real<lower = 0> omega_ec50 = omega_pd[3];

  real<lower = 0> omega_sq_e0 = square(omega_e0);
  real<lower = 0> omega_sq_ke0 = square(omega_ke0);
  real<lower = 0> omega_sq_ec50 = square(omega_ec50);

  real cor_e0_ke0;
  real cor_e0_ec50;
  real cor_ke0_ec50;
  
  real omega_e0_ke0;
  real omega_e0_ec50;
  real omega_ke0_ec50;

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
    cor_cl_ka = R[1, 3];
    cor_vc_ka = R[2, 3];

    omega_cl_vc = Omega[1, 2];
    omega_cl_ka = Omega[1, 3];
    omega_vc_ka = Omega[2, 3];
    
    cor_e0_ke0 = R_pd[1, 2];
    cor_e0_ec50 = R_pd[1, 3];
    cor_ke0_ec50 = R_pd[2, 3];

    omega_e0_ke0 = Omega_pd[1, 2];
    omega_e0_ec50 = Omega_pd[1, 3];
    omega_ke0_ec50 = Omega_pd[2, 3];
  
  }

  if(no_gq_predictions == 0){
    
    vector[n_total] dv_pred;
    matrix[n_total, n_cmt + n_cmt_pd] x_pred;
    vector[n_total] dv_ipred;
    matrix[n_total, n_cmt + n_cmt_pd] x_ipred;

    for(j in 1:n_subjects){
        
      matrix[n_cmt + n_cmt_pd, n_cmt + n_cmt_pd] K = 
                            rep_matrix(0, n_cmt + n_cmt_pd, n_cmt + n_cmt_pd);
                            
      matrix[n_cmt + n_cmt_pd, n_cmt + n_cmt_pd] K_tv = 
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
                         
      K_tv[1, 1] = -TVKA;
      K_tv[2, 1] = TVKA;
      K_tv[2, 2] = -TVCL/TVVC;
      K_tv[3, 2] = TVKE0;
      K_tv[3, 3] = -TVKE0;

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
      
      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){
          dv_ipred[k] = x_ipred[k, 2] / VC[j];
          dv_pred[k] = x_pred[k, 2] / TVVC;
        }else if(cmt[k] == 3){
          real conc_effcmt = x_ipred[k, 3]/VC[j];
          real inh = (EMAX[j]*pow(conc_effcmt, HILL[j]))/(pow(EC50[j], HILL[j]) + 
                                                    pow(conc_effcmt, HILL[j]));
          dv_ipred[k] = E0[j]*(1 - inh);
          
          real tv_conc_effcmt = x_pred[k, 3]/TVVC;
          real tvinh = (TVEMAX*pow(tv_conc_effcmt, TVHILL))/(pow(TVEC50, TVHILL) + 
                                                    pow(tv_conc_effcmt, TVHILL));
          dv_pred[k] = TVE0*(1 - tvinh);
        }
      }
    }
    
    pred = dv_pred[i_obs];
    ipred = dv_ipred[i_obs];
      
    res = log(dv_obs) - log(pred);
    ires = log(dv_obs) - log(ipred);

    for(i in 1:n_obs){
      
      if(cmt[i_obs[i]] == 2 || cmt[i_obs[i]] == 3){
        real log_ipred_tmp = log(ipred[i]);
        real sigma_tmp = cmt[i_obs[i]] == 2 ? sigma : sigma_pd;
      
        dv_ppc[i] = lognormal_rng(log_ipred_tmp, sigma);
      
        if(bloq_obs[i] == 1){
          log_lik[i] = lognormal_lcdf(lloq_obs[i] | log_ipred_tmp, sigma_tmp);
        }else{
          log_lik[i] = lognormal_lpdf(dv_obs[i] | log_ipred_tmp, sigma_tmp);
        }
      
        wres[i] = res[i]/sigma_tmp;
        iwres[i] = ires[i]/sigma_tmp;
      }
    }
  }
}

