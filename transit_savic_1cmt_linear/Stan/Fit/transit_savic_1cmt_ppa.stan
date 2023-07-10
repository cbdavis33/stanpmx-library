// One-compartment PK Model with Transit Compartment absorption
// IIV on CL, VC, KA, NTR, MTT (full covariance matrix)
// n_transit is a positive real number, not fixed, and not necessarily an integer
// General ODE solution using pure Stan code
// proportional plus additive error - DV = IPRED*(1 + eps_p) + eps_a
// Implements threading for within-chain parallelization
// Deals with BLOQ values by the "CDF trick" (M4)
// The 0th transit compartment is where the dosing happens (cmt = 1). This could
//   also be called the Depot
// Drug is absorbed from the depot through the transit compartments into the 
//   absorption compartment and into the central compartment
//   (depot -> tr1 -> ... -> trn -> absorption -> central)
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
  
  vector transit_1cmt_ode(real t, 
                          vector y,
                          real ke, real ka, real n_tr, real f, 
                          real ktr,
                          int n_dose_subj,            // # of doses for this subject
                          array[] real dosetime_subj, // dosetimes for this subject
                          array[] real doseamt_subj){ // dose amounts for this subject
    
    // real ktr = (n_tr + 1)/mtt;
    real k_inpt = f*pow(ktr, n_tr + 1)/exp(lgamma(n_tr + 1));
    
    real inpt = 0;
    array[n_dose_subj] real ipt = rep_array(0.0, n_dose_subj);
    
    for(i in 1:n_dose_subj){
      if(t >= dosetime_subj[i]){
        real delta_t = t - dosetime_subj[i];
        ipt[i] = doseamt_subj[i]*pow(delta_t, n_tr)*exp(-ktr*delta_t);
      }
    }
    
    inpt = sum(ipt);
    
    vector[2] dydt;
    
    dydt[1] = k_inpt*inpt - ka*y[1];
    dydt[2] = ka*y[1] - ke*y[2];
    
    return dydt;
    
  }
  
  real partial_sum_lpmf(array[] int seq_subj, int start, int end,
                        vector dv_obs, array[] int dv_obs_id, array[] int i_obs,
                        array[] real amt, array[] int cmt, array[] int evid, 
                        array[] real time, array[] real rate, array[] real ii, 
                        array[] int addl, array[] int ss,
                        array[] int subj_start, array[] int subj_end,
                        array[] int subj_start_dose, array[] int subj_end_dose,
                        vector KE, vector VC, vector KA, vector NTR, vector KTR,
                        real sigma_sq_p, real sigma_sq_a, real sigma_p_a,
                        vector lloq, array[] int bloq,
                        int n_random, int n_subjects, int n_total,
                        array[] real dosetime, array[] real doseamt,
                        array[] real bioav, int n_cmt, vector y0, int solver){
                           
    real ptarget = 0;
    
    int N = end - start + 1;        // number of subjects in this slice  
    
    vector[n_total] dv_ipred;       ;
    array[n_total] vector[n_cmt] x_ipred;

    int n_obs_slice = num_between(subj_start[start], subj_end[end], i_obs);
    array[n_obs_slice] int i_obs_slice = find_between(subj_start[start], 
                                                      subj_end[end], i_obs);
                                                
    vector[n_obs_slice] dv_obs_slice = find_between_vec(start, end, 
                                                        dv_obs_id, dv_obs);
    
    vector[n_obs_slice] ipred_slice;
    
    vector[n_obs_slice] lloq_slice = lloq[i_obs_slice];
    array[n_obs_slice] int bloq_slice = bloq[i_obs_slice];
    
    for(n in 1:N){            // loop over subjects in this slice
    
      int j = n + start - 1; // j is the ID of the current subject
      
      real t0 = time[subj_start[j]];
      int n_dose_subj = subj_end_dose[j] - subj_start_dose[j] + 1;
      
      array[n_dose_subj] real dosetime_subj = 
        dosetime[subj_start_dose[j]:subj_end_dose[j]];
      array[n_dose_subj] real doseamt_subj = 
        doseamt[subj_start_dose[j]:subj_end_dose[j]];
        
      x_ipred[subj_start[j],] = y0;

      if(solver == 1){
        x_ipred[(subj_start[j] + 1):subj_end[j],] = 
          ode_rk45(transit_1cmt_ode, y0, t0, 
                   time[(subj_start[j] + 1):subj_end[j]], 
                   KE[j], KA[j], NTR[j], bioav[j], KTR[j],
                   n_dose_subj, dosetime_subj, doseamt_subj);
      }else if(solver == 2){
        x_ipred[(subj_start[j] + 1):subj_end[j],] = 
          ode_bdf(transit_1cmt_ode, y0, t0, 
                  time[(subj_start[j] + 1):subj_end[j]], 
                  KE[j], KA[j], NTR[j], bioav[j], KTR[j],
                  n_dose_subj, dosetime_subj, doseamt_subj);
      }else if(solver == 3){
        x_ipred[(subj_start[j] + 1):subj_end[j],] = 
          ode_adams(transit_1cmt_ode, y0, t0, 
                   time[(subj_start[j] + 1):subj_end[j]], 
                   KE[j], KA[j], NTR[j], bioav[j], KTR[j],
                   n_dose_subj, dosetime_subj, doseamt_subj);
      }else{
        x_ipred[(subj_start[j] + 1):subj_end[j],] = 
          ode_ckrk(transit_1cmt_ode, y0, t0, 
                   time[(subj_start[j] + 1):subj_end[j]], 
                   KE[j], KA[j], NTR[j], bioav[j], KTR[j],
                   n_dose_subj, dosetime_subj, doseamt_subj);
      }
      
      for(k in subj_start[j]:subj_end[j]){
        dv_ipred[k] = fmax(0.000000001, x_ipred[k, 2] / VC[j]);
        // dv_ipred[k] = x_ipred[k, 2] / VC[j];
        
      }
    
    }
  
    ipred_slice = dv_ipred[i_obs_slice];
    
    for(i in 1:n_obs_slice){
      
      real ipred_tmp = ipred_slice[i];
      real sigma_tmp = sqrt(square(ipred_tmp) * sigma_sq_p + sigma_sq_a + 
                            2*ipred_tmp*sigma_p_a);
      
      if(bloq_slice[i] == 1){
        ptarget += log_diff_exp(normal_lcdf(lloq_slice[i] | ipred_tmp, 
                                                            sigma_tmp),
                                normal_lcdf(0.0 | ipred_tmp, sigma_tmp)) -
                   normal_lccdf(0.0 | ipred_tmp, sigma_tmp); 
      }else{
        ptarget += normal_lpdf(dv_obs_slice[i] | ipred_tmp, sigma_tmp) -
                   normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
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
  
  real<lower = 0> location_tvcl;   // Prior Location parameter for CL
  real<lower = 0> location_tvvc;   // Prior Location parameter for VC
  real<lower = 0> location_tvka;   // Prior Location parameter for KA
  real<lower = 0> location_tvntr;  // Prior Location parameter for NTR
  real<lower = 0> location_tvmtt;  // Prior Location parameter for MTT
  
  real<lower = 0> scale_tvcl;      // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;      // Prior Scale parameter for VC
  real<lower = 0> scale_tvka;      // Prior Scale parameter for KA
  real<lower = 0> scale_tvntr;     // Prior Scale parameter for NTR
  real<lower = 0> scale_tvmtt;     // Prior Scale parameter for MTT
  
  real<lower = 0> scale_omega_cl;  // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc;  // Prior scale parameter for omega_vc
  real<lower = 0> scale_omega_ka;  // Prior scale parameter for omega_ka
  real<lower = 0> scale_omega_ntr; // Prior scale parameter for omega_ntr
  real<lower = 0> scale_omega_mtt; // Prior scale parameter for omega_mtt
  
  real<lower = 0> lkj_df_omega;   // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma_p;  // Prior Scale parameter for proportional error
  real<lower = 0> scale_sigma_a;  // Prior Scale parameter for additive error
  
  real<lower = 0> lkj_df_sigma;   // Prior degrees of freedom for sigma cor mat
  
  int n_dose;
  array[n_dose] real dosetime;
  array[n_dose] real doseamt;
  
  array[n_subjects] int subj_start_dose;
  array[n_subjects] int subj_end_dose;
  
  int<lower = 0, upper = 1> prior_only;
  
  int<lower = 1, upper = 4> solver; // 1 = rk45, 2 = bdf, 3 = adams, 4 = ckrk
 
}
transformed data{ 
  
  int grainsize = 1;
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 5;           // Number of random effects
  int n_cmt = 2;              // Absorption, central
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc,
                                      scale_omega_ka, scale_omega_ntr,
                                      scale_omega_mtt}; 
                                      
  array[2] real scale_sigma = {scale_sigma_p, scale_sigma_a}; 
  
  array[n_subjects] int seq_subj = sequence(1, n_subjects); // reduce_sum over subjects
  
  array[n_subjects] real bioav = rep_array(1.0, n_subjects); // Hardcoding, but could be data or a parameter in another situation

  vector[n_cmt] y0 = rep_vector(0.0, n_cmt);
  
}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = TVCL/TVVC> TVKA;
  real<lower = 0> TVNTR; 
  real<lower = 0> TVMTT; 
  
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
  vector[n_subjects] eta_ntr;
  vector[n_subjects] eta_mtt;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  vector[n_subjects] NTR;
  vector[n_subjects] MTT;
  vector[n_subjects] KE;
  vector[n_subjects] KTR;
  
  real<lower = 0> sigma_p = sigma[1];
  real<lower = 0> sigma_a = sigma[2];
  
  real<lower = 0> sigma_sq_p = square(sigma_p);
  real<lower = 0> sigma_sq_a = square(sigma_a);
  
  real cor_p_a;
  real sigma_p_a;

  {
  
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA, 
                                                         TVNTR, TVMTT});

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
    
    eta_cl = col(eta, 1);
    eta_vc = col(eta, 2);
    eta_ka = col(eta, 3);
    eta_ntr = col(eta, 4);
    eta_mtt = col(eta, 5);
    CL = col(theta, 1);
    VC = col(theta, 2);
    KA = col(theta, 3);
    NTR = col(theta, 4);
    MTT = col(theta, 5);
    KE = CL ./ VC;
    KTR = (NTR + 1) ./ MTT;
    
    cor_p_a = R_Sigma[1, 2];
    sigma_p_a = Sigma[1, 2];
  
  }
  
}
model{ 
  
  // Priors
  TVCL ~ lognormal(log(location_tvcl), scale_tvcl);
  TVVC ~ lognormal(log(location_tvvc), scale_tvvc);
  TVKA ~ lognormal(log(location_tvka), scale_tvka) T[TVCL/TVVC, ];
  TVNTR ~ lognormal(log(location_tvntr), scale_tvntr); 
  TVMTT ~ lognormal(log(location_tvmtt), scale_tvmtt);

  omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma ~ normal(0, scale_sigma);
  L_Sigma ~ lkj_corr_cholesky(lkj_df_sigma);
  
  to_vector(Z) ~ std_normal();
  
  // Likelihood
  if(prior_only == 0){
    target += reduce_sum(partial_sum_lupmf, seq_subj, grainsize,
                         dv_obs, dv_obs_id, i_obs,
                         amt, cmt, evid, time, 
                         rate, ii, addl, ss, 
                         subj_start, subj_end, subj_start_dose, subj_end_dose,
                         KE, VC, KA, NTR, KTR,
                         sigma_sq_p, sigma_sq_a, sigma_p_a, 
                         lloq, bloq,
                         n_random, n_subjects, n_total,
                         dosetime, doseamt,
                         bioav, n_cmt, y0, solver);
  } 
}
generated quantities{

  real<lower = 0> TVKE = TVCL/TVVC;
  real<lower = 0> TVKTR = (TVNTR + 1)/TVMTT;

  real<lower = 0> omega_cl = omega[1];
  real<lower = 0> omega_vc = omega[2];
  real<lower = 0> omega_ka = omega[3];
  real<lower = 0> omega_ntr = omega[4];
  real<lower = 0> omega_mtt = omega[5];

  real<lower = 0> omega_sq_cl = square(omega_cl);
  real<lower = 0> omega_sq_vc = square(omega_vc);
  real<lower = 0> omega_sq_ka = square(omega_ka);
  real<lower = 0> omega_sq_ntr = square(omega_ntr);
  real<lower = 0> omega_sq_mtt = square(omega_mtt);

  real cor_cl_vc;
  real cor_cl_ka;
  real cor_cl_ntr;
  real cor_cl_mtt;
  real cor_vc_ka;
  real cor_vc_ntr;
  real cor_vc_mtt;
  real cor_ka_ntr;
  real cor_ka_mtt;
  real cor_ntr_mtt;

  real omega_cl_vc;
  real omega_cl_ka;
  real omega_cl_ntr;
  real omega_cl_mtt;
  real omega_vc_ka;
  real omega_vc_ntr;
  real omega_vc_mtt;
  real omega_ka_ntr;
  real omega_ka_mtt;
  real omega_ntr_mtt;

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

    vector[n_total] dv_ipred;
    array[n_total] vector[n_cmt] x_ipred;
    vector[n_total] dv_pred;
    array[n_total] vector[n_cmt] x_pred;

    cor_cl_vc = R[1, 2];
    cor_cl_ka = R[1, 3];
    cor_cl_ntr = R[1, 4];
    cor_cl_mtt = R[1, 5];
    cor_vc_ka = R[2, 3];
    cor_vc_ntr = R[2, 4];
    cor_vc_mtt = R[2, 5];
    cor_ka_ntr = R[3, 4];
    cor_ka_mtt = R[3, 5];
    cor_ntr_mtt = R[4, 5];

    omega_cl_vc = Omega[1, 2];
    omega_cl_ka = Omega[1, 3];
    omega_cl_ntr = Omega[1, 4];
    omega_cl_mtt = Omega[1, 5];
    omega_vc_ka = Omega[2, 3];
    omega_vc_ntr = Omega[2, 4];
    omega_vc_mtt = Omega[2, 5];
    omega_ka_ntr = Omega[3, 4];
    omega_ka_mtt = Omega[3, 5];
    omega_ntr_mtt = Omega[4, 5];

    for(j in 1:n_subjects){

      real t0 = time[subj_start[j]];
      int n_dose_subj = subj_end_dose[j] - subj_start_dose[j] + 1;

      array[n_dose_subj] real dosetime_subj =
        dosetime[subj_start_dose[j]:subj_end_dose[j]];
      array[n_dose_subj] real doseamt_subj =
        doseamt[subj_start_dose[j]:subj_end_dose[j]];

      x_ipred[subj_start[j],] = y0;
      x_pred[subj_start[j],] = y0;

      if(solver == 1){
        x_ipred[(subj_start[j] + 1):subj_end[j],] =
          ode_rk45(transit_1cmt_ode, y0, t0,
                   time[(subj_start[j] + 1):subj_end[j]],
                   KE[j], KA[j], NTR[j], bioav[j], KTR[j],
                   n_dose_subj, dosetime_subj, doseamt_subj);

        x_pred[(subj_start[j] + 1):subj_end[j],] =
          ode_rk45(transit_1cmt_ode, y0, t0,
                   time[(subj_start[j] + 1):subj_end[j]],
                   TVKE, TVKA, TVNTR, bioav[j], TVKTR,
                   n_dose_subj, dosetime_subj, doseamt_subj);
      }else if(solver == 2){
        x_ipred[(subj_start[j] + 1):subj_end[j],] =
          ode_bdf(transit_1cmt_ode, y0, t0,
                   time[(subj_start[j] + 1):subj_end[j]],
                   KE[j], KA[j], NTR[j], bioav[j], KTR[j],
                   n_dose_subj, dosetime_subj, doseamt_subj);

        x_pred[(subj_start[j] + 1):subj_end[j],] =
          ode_bdf(transit_1cmt_ode, y0, t0,
                  time[(subj_start[j] + 1):subj_end[j]],
                  TVKE, TVKA, TVNTR, bioav[j], TVKTR,
                  n_dose_subj, dosetime_subj, doseamt_subj);
      }else if(solver == 3){
        x_ipred[(subj_start[j] + 1):subj_end[j],] =
          ode_adams(transit_1cmt_ode, y0, t0,
                    time[(subj_start[j] + 1):subj_end[j]],
                    KE[j], KA[j], NTR[j], bioav[j], KTR[j],
                    n_dose_subj, dosetime_subj, doseamt_subj);

        x_pred[(subj_start[j] + 1):subj_end[j],] =
          ode_adams(transit_1cmt_ode, y0, t0,
                    time[(subj_start[j] + 1):subj_end[j]],
                    TVKE, TVKA, TVNTR, bioav[j], TVKTR,
                    n_dose_subj, dosetime_subj, doseamt_subj);
      }else{
        x_ipred[(subj_start[j] + 1):subj_end[j],] =
          ode_ckrk(transit_1cmt_ode, y0, t0,
                   time[(subj_start[j] + 1):subj_end[j]],
                   KE[j], KA[j], NTR[j], bioav[j], KTR[j],
                   n_dose_subj, dosetime_subj, doseamt_subj);

        x_pred[(subj_start[j] + 1):subj_end[j],] =
          ode_ckrk(transit_1cmt_ode, y0, t0,
                   time[(subj_start[j] + 1):subj_end[j]],
                   TVKE, TVKA, TVNTR, bioav[j], TVKTR,
                   n_dose_subj, dosetime_subj, doseamt_subj);
      }

      for(k in subj_start[j]:subj_end[j]){
        dv_ipred[k] = fmax(0.000000001, x_ipred[k, 2] / VC[j]);
        dv_pred[k] = fmax(0.000000001, x_pred[k, 2] / TVVC);
      }

    }

    pred = dv_pred[i_obs];
    ipred = dv_ipred[i_obs];

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





