// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model with linear and Michaelis-Menten elimination
// IIV on CL, VC, VMAX, KM, KA (full covariance matrix)
// exponential error - DV = IPRED*exp(eps)
// General ODE solution (rk45) using Torsten 
// Implements threading for within-chain parallelization 
// Deals with BLOQ values by the M3 method (M3 and M4 are equivalent with this
//   error model)

functions{
  
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
  
  vector find_between(int lb, int ub, array[] int idx, vector y) {
    
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
  
  vector depot_1cmt_linear_plus_mm_ode(real t, vector y, array[] real params, 
                                       array[] real x_r, array[] int x_i){
  
    real cl = params[1];
    real vc = params[2];
    real vmax = params[3];
    real km = params[4];
    real ka = params[5];
    real ke = cl/vc;
    
    real conc = y[2]/vc;
    
    vector[2] dydt;

    dydt[1] = -ka*y[1];
    dydt[2] = ka*y[1] - ke*y[2] - vmax*conc/(km + conc);
    
    return dydt;
    
  }
  
  real partial_sum_lpmf(array[] int seq_subj, int start, int end,
                        vector dv_obs, array[] int dv_obs_id, array[] int i_obs,
                        array[] real amt, array[] int cmt, array[] int evid, 
                        array[] real time, array[] real rate, array[] real ii, 
                        array[] int addl, array[] int ss,
                        array[] int subj_start, array[] int subj_end, 
                        vector CL, vector VC, vector VMAX, vector KM, vector KA, 
                        real sigma, 
                        vector lloq, array[] int bloq,
                        int n_random, int n_subjects, int n_total,
                        array[] real bioav, array[] real tlag, int n_cmt){
                           
    real ptarget = 0;
                              
    int N = end - start + 1;    // number of subjects in this slice  
    vector[n_total] dv_ipred;   
    matrix[n_total, 2] x_ipred;
  
    int n_obs_slice = num_between(subj_start[start], subj_end[end], i_obs);
    array[n_obs_slice] int i_obs_slice = find_between(subj_start[start], 
                                                      subj_end[end], i_obs);
                                                
    vector[n_obs_slice] dv_obs_slice = find_between(start, end, 
                                                    dv_obs_id, dv_obs);
    
    vector[n_obs_slice] ipred_slice;
    
    vector[n_obs_slice] lloq_slice = lloq[i_obs_slice];
    array[n_obs_slice] int bloq_slice = bloq[i_obs_slice];
    
    
    for(n in 1:N){            // loop over subjects in this slice
    
      int j = n + start - 1; // j is the ID of the current subject
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_1cmt_linear_plus_mm_ode,
                       n_cmt,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {CL[j], VC[j], VMAX[j], KM[j], KA[j]},
                       bioav, tlag)';
                       
      
      dv_ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
    
    }
  
    ipred_slice = dv_ipred[i_obs_slice];
    
    for(i in 1:n_obs_slice){
      if(bloq_slice[i] == 1){
        ptarget += lognormal_lcdf(lloq_slice[i] | log(ipred_slice[i]), sigma);
      }else{
        ptarget += lognormal_lpdf(dv_obs_slice[i] | log(ipred_slice[i]), sigma);
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
  
  real<lower = 0> location_tvcl;    // Prior Location parameter for CL
  real<lower = 0> location_tvvc;    // Prior Location parameter for VC
  real<lower = 0> location_tvvmax;  // Prior Location parameter for VMAX
  real<lower = 0> location_tvkm;    // Prior Location parameter for KM
  real<lower = 0> location_tvka;    // Prior Location parameter for KA
  
  real<lower = 0> scale_tvcl;       // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;       // Prior Scale parameter for VC
  real<lower = 0> scale_tvvmax;     // Prior Scale parameter for VMAX
  real<lower = 0> scale_tvkm;       // Prior Scale parameter for KM
  real<lower = 0> scale_tvka;       // Prior Scale parameter for KA
  
  real<lower = 0> scale_omega_cl;   // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc;   // Prior scale parameter for omega_vc
  real<lower = 0> scale_omega_vmax; // Prior scale parameter for omega_vmax
  real<lower = 0> scale_omega_km;   // Prior scale parameter for omega_km
  real<lower = 0> scale_omega_ka;   // Prior scale parameter for omega_ka
  
  real<lower = 0> lkj_df_omega;   // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma;    // Prior Scale parameter for lognormal error
  
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
  
  int n_random = 5;                    // Number of random effects
  int n_cmt = 2;                       // Number of states in the ODEs
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc, 
                                      scale_omega_vmax, scale_omega_km, 
                                      scale_omega_ka};
  
  array[n_subjects] int seq_subj = linspaced_int_array(n_subjects, 1, n_subjects); // reduce_sum over subjects
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
}
parameters{ 
  
  real<lower = 0> TVCL;
  real<lower = 0> TVVC;
  real<lower = 0> TVVMAX;       
  real<lower = 0> TVKM;
  real<lower = 0> TVKA;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
transformed parameters{
  
  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_vc;
  vector[n_subjects] eta_vmax;
  vector[n_subjects] eta_km;
  vector[n_subjects] eta_ka;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] VMAX;
  vector[n_subjects] KM;
  vector[n_subjects] KA;

  {
  
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVVMAX, TVKM, TVKA});

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
    
    eta_cl = col(eta, 1);
    eta_vc = col(eta, 2);
    eta_vmax = col(eta, 3);
    eta_km = col(eta, 4);
    eta_ka = col(eta, 5);
    CL = col(theta, 1);
    VC = col(theta, 2);
    VMAX = col(theta, 3);
    KM = col(theta, 4);
    KA = col(theta, 5);
  
  }
  
}
model{ 
  
  // Priors
  TVCL ~ lognormal(log(location_tvcl), scale_tvcl);
  TVVC ~ lognormal(log(location_tvvc), scale_tvvc);
  TVVMAX ~ lognormal(log(location_tvvmax), scale_tvvmax);
  TVKM ~ lognormal(log(location_tvkm), scale_tvkm);
  TVKA ~ lognormal(log(location_tvka), scale_tvka);

  omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma ~ normal(0, scale_sigma);
  
  to_vector(Z) ~ std_normal();
  
  // Likelihood
  if(prior_only == 0){
    target += reduce_sum(partial_sum_lupmf, seq_subj, grainsize,
                         dv_obs, dv_obs_id, i_obs,
                         amt, cmt, evid, time, 
                         rate, ii, addl, ss, subj_start, subj_end, 
                         CL, VC, VMAX, KM, KA,
                         sigma,
                         lloq, bloq,
                         n_random, n_subjects, n_total,
                         bioav, tlag, n_cmt);
  }
}
generated quantities{
  
  real<lower = 0> sigma_sq = square(sigma);

  real<lower = 0> omega_cl = omega[1];
  real<lower = 0> omega_vc = omega[2];
  real<lower = 0> omega_vmax = omega[3];
  real<lower = 0> omega_km = omega[4];
  real<lower = 0> omega_ka = omega[5];

  real<lower = 0> omega_sq_cl = square(omega_cl);
  real<lower = 0> omega_sq_vc = square(omega_vc);
  real<lower = 0> omega_sq_vmax = square(omega_vmax);
  real<lower = 0> omega_sq_km = square(omega_km);
  real<lower = 0> omega_sq_ka = square(omega_ka);

  real cor_cl_vc;
  real cor_cl_vmax;
  real cor_cl_km;
  real cor_cl_ka;
  real cor_vc_vmax;
  real cor_vc_km;
  real cor_vc_ka;
  real cor_vmax_km;
  real cor_vmax_ka;
  real cor_km_ka;
  
  real omega_cl_vc;
  real omega_cl_vmax;
  real omega_cl_km;
  real omega_cl_ka;
  real omega_vc_vmax;
  real omega_vc_km;
  real omega_vc_ka;
  real omega_vmax_km;
  real omega_vmax_ka;
  real omega_km_ka;

  vector[no_gq_predictions ? 0 : n_obs] pred;
  vector[no_gq_predictions ? 0 : n_obs] epred_stan;
  vector[no_gq_predictions ? 0 : n_obs] ipred;
  vector[no_gq_predictions ? 0 : n_obs] epred;
  vector[no_gq_predictions ? 0 : n_obs] dv_ppc;
  vector[no_gq_predictions ? 0 : n_obs] log_lik;
  vector[no_gq_predictions ? 0 : n_obs] iwres;
 
  {

    matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
    matrix[n_random, n_random] Omega = quad_form_diag(R, omega);

    cor_cl_vc = R[1, 2];
    cor_cl_vmax = R[1, 3];
    cor_cl_km = R[1, 4];
    cor_cl_ka = R[1, 5];
    cor_vc_vmax = R[2, 3];
    cor_vc_km = R[2, 4];
    cor_vc_ka = R[2, 5];
    cor_vmax_km = R[3, 4];
    cor_vmax_ka = R[3, 5];
    cor_km_ka = R[4, 5];
    
    omega_cl_vc = Omega[1, 2];
    omega_cl_vmax = Omega[1, 3];
    omega_cl_km = Omega[1, 4];
    omega_cl_ka = Omega[1, 5];
    omega_vc_vmax = Omega[1, 2];
    omega_vc_km = Omega[1, 3];
    omega_vc_ka = Omega[1, 4];
    omega_vmax_km = Omega[2, 3];
    omega_vmax_ka = Omega[2, 4];
    omega_km_ka = Omega[3, 4];
    
  }

  if(no_gq_predictions == 0){
    
    vector[n_subjects] CL_new;
    vector[n_subjects] VC_new;
    vector[n_subjects] VMAX_new;
    vector[n_subjects] KM_new;
    vector[n_subjects] KA_new;
    
    vector[n_total] dv_pred;
    matrix[n_total, n_cmt] x_pred;
    vector[n_total] dv_epred;
    matrix[n_total, n_cmt] x_epred;
    vector[n_total] dv_ipred;
    matrix[n_total, n_cmt] x_ipred;
    
    matrix[n_subjects, n_random] eta_new;
    matrix[n_subjects, n_random] theta_new;
    
    for(i in 1:n_subjects){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(to_row_vector({TVCL, TVVC, TVVMAX, TVKM, TVKA}), n_subjects) .* 
                                                                  exp(eta_new));
    
    for(j in 1:n_subjects){
      
      row_vector[n_random] theta_j_new = theta_new[j]; // access the parameters for subject j's epred
      
      real cl_p = TVCL;
      real vc_p = TVVC;
      real vmax_p = TVVMAX;
      real km_p = TVKM;
      real ka_p = TVKA;
      
      CL_new[j] = theta_j_new[1];
      VC_new[j] = theta_j_new[2];
      VMAX_new[j] = theta_j_new[3];
      KM_new[j] = theta_j_new[4];
      KA_new[j] = theta_j_new[5];
      
      x_ipred[subj_start[j]:subj_end[j],] =
            pmx_solve_rk45(depot_1cmt_linear_plus_mm_ode,
                           n_cmt,
                           time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           {CL[j], VC[j], VMAX[j], KM[j], KA[j]},
                           bioav, tlag)';

      x_epred[subj_start[j]:subj_end[j],] =
            pmx_solve_rk45(depot_1cmt_linear_plus_mm_ode,
                           n_cmt,
                           time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           {CL_new[j], VC_new[j], VMAX_new[j], KM_new[j], KA_new[j]},
                           bioav, tlag)';
                           
      x_pred[subj_start[j]:subj_end[j],] =
            pmx_solve_rk45(depot_1cmt_linear_plus_mm_ode,
                           n_cmt,
                           time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           {cl_p, vc_p, vmax_p, km_p, ka_p},
                           bioav, tlag)';
                           
      dv_ipred[subj_start[j]:subj_end[j]] =
        x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
        
      dv_epred[subj_start[j]:subj_end[j]] =
        x_epred[subj_start[j]:subj_end[j], 2] ./ VC_new[j];
      
      dv_pred[subj_start[j]:subj_end[j]] =
        x_pred[subj_start[j]:subj_end[j], 2] ./ vc_p;
    }

    pred = dv_pred[i_obs];
    epred_stan = dv_epred[i_obs];
    ipred = dv_ipred[i_obs];

    for(i in 1:n_obs){
      real log_ipred_tmp = log(ipred[i]);
      real log_epred_tmp = log(epred_stan[i]);
      dv_ppc[i] = lognormal_rng(log_ipred_tmp, sigma);
      epred[i] = lognormal_rng(log_epred_tmp, sigma);
      if(bloq_obs[i] == 1){
        log_lik[i] = lognormal_lcdf(lloq_obs[i] | log_ipred_tmp, sigma);
      }else{
        log_lik[i] = lognormal_lpdf(dv_obs[i] | log_ipred_tmp, sigma);
      }
    }
    iwres = (log(dv_obs) - log(ipred))/sigma;
  }
}
