// Two-compartment PK Model with Transit Compartment absorption
// IIV on CL, VC, Q, VP, KA, NTR, MTT (full covariance matrix)
// n_transit is a positive real number, not fixed, and not necessarily an integer
// General ODE solution using pure Stan code
// proportional error - DV = IPRED*(1 + eps_p)
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
  
  vector transit_2cmt_ode(real t, 
                          vector y,
                          real cl, real vc, real q, real vp, real ka, 
                          real n_tr, real f, real ktr, 
                          int n_dose_subj,            // # of doses for this subject
                          array[] real dosetime_subj, // dosetimes for this subject
                          array[] real doseamt_subj){ // dose amounts for this subject
    
    // real ktr = (n_tr + 1)/mtt;
    real k_inpt = f*pow(ktr, n_tr + 1)/exp(lgamma(n_tr + 1));
    
    real inpt = 0;
    array[n_dose_subj] real ipt = rep_array(0.0, n_dose_subj);
    
    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;
    
    vector[3] dydt;
    
    for(i in 1:n_dose_subj){
      if(t >= dosetime_subj[i]){
        real delta_t = t - dosetime_subj[i];
        ipt[i] = doseamt_subj[i]*pow(delta_t, n_tr)*exp(-ktr*delta_t);
      }
    }
    
    inpt = sum(ipt);
    
    dydt[1] = k_inpt*inpt - ka*y[1];
    dydt[2] = ka*y[1] - (ke + k_cp)*y[2] + k_pc*y[3];
    dydt[3] = k_cp*y[2] - k_pc*y[3];
    
    return dydt;
    
  }
  
  real partial_sum_lpmf(array[] int seq_subj, int start, int end,
                        vector dv_obs, array[] int dv_obs_id, array[] int i_obs,
                        array[] real amt, array[] int cmt, array[] int evid, 
                        array[] real time, array[] real rate, array[] real ii, 
                        array[] int addl, array[] int ss,
                        array[] int subj_start, array[] int subj_end,
                        array[] int subj_start_dose, array[] int subj_end_dose,
                        vector CL, vector VC, vector Q, vector VP, vector KA, 
                        vector NTR, vector KTR,
                        real sigma_p,
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
          ode_rk45(transit_2cmt_ode, y0, t0, 
                   time[(subj_start[j] + 1):subj_end[j]], 
                   CL[j], VC[j], Q[j], VP[j], KA[j], NTR[j], bioav[j], KTR[j],
                   n_dose_subj, dosetime_subj, doseamt_subj);
      }else if(solver == 2){
        x_ipred[(subj_start[j] + 1):subj_end[j],] = 
          ode_bdf(transit_2cmt_ode, y0, t0, 
                  time[(subj_start[j] + 1):subj_end[j]], 
                  CL[j], VC[j], Q[j], VP[j], KA[j], NTR[j], bioav[j], KTR[j],
                  n_dose_subj, dosetime_subj, doseamt_subj);
      }else if(solver == 3){
        x_ipred[(subj_start[j] + 1):subj_end[j],] = 
          ode_adams(transit_2cmt_ode, y0, t0, 
                   time[(subj_start[j] + 1):subj_end[j]], 
                   CL[j], VC[j], Q[j], VP[j], KA[j], NTR[j], bioav[j], KTR[j],
                   n_dose_subj, dosetime_subj, doseamt_subj);
      }else{
        x_ipred[(subj_start[j] + 1):subj_end[j],] = 
          ode_ckrk(transit_2cmt_ode, y0, t0, 
                   time[(subj_start[j] + 1):subj_end[j]], 
                   CL[j], VC[j], Q[j], VP[j], KA[j], NTR[j], bioav[j], KTR[j],
                   n_dose_subj, dosetime_subj, doseamt_subj);
      }
      
      for(k in subj_start[j]:subj_end[j]){
        dv_ipred[k] = fmax(1e-14, x_ipred[k, 2] / VC[j]);
      }
    
    }
  
    ipred_slice = dv_ipred[i_obs_slice];
    
    for(i in 1:n_obs_slice){
      real sigma_tmp = ipred_slice[i]*sigma_p;
      if(bloq_slice[i] == 1){
        ptarget += log_diff_exp(normal_lcdf(lloq_slice[i] | ipred_slice[i], 
                                                            sigma_tmp),
                                normal_lcdf(0.0 | ipred_slice[i], sigma_tmp)) -
                   normal_lccdf(0.0 | ipred_slice[i], sigma_tmp); 
      }else{
        ptarget += normal_lpdf(dv_obs_slice[i] | ipred_slice[i], sigma_tmp) -
                   normal_lccdf(0.0 | ipred_slice[i], sigma_tmp);
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
  real<lower = 0> location_tvq;    // Prior Location parameter for Q
  real<lower = 0> location_tvvp;   // Prior Location parameter for VP
  real<lower = 0> location_tvka;   // Prior Location parameter for KA
  real<lower = 0> location_tvntr;  // Prior Location parameter for NTR
  real<lower = 0> location_tvmtt;  // Prior Location parameter for MTT
  
  real<lower = 0> scale_tvcl;      // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;      // Prior Scale parameter for VC
  real<lower = 0> scale_tvq;       // Prior Scale parameter for Q
  real<lower = 0> scale_tvvp;      // Prior Scale parameter for VP
  real<lower = 0> scale_tvka;      // Prior Scale parameter for KA
  real<lower = 0> scale_tvntr;     // Prior Scale parameter for NTR
  real<lower = 0> scale_tvmtt;     // Prior Scale parameter for MTT
  
  real<lower = 0> scale_omega_cl;  // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc;  // Prior scale parameter for omega_vc
  real<lower = 0> scale_omega_q;   // Prior scale parameter for omega_q
  real<lower = 0> scale_omega_vp;  // Prior scale parameter for omega_vp
  real<lower = 0> scale_omega_ka;  // Prior scale parameter for omega_ka
  real<lower = 0> scale_omega_ntr; // Prior scale parameter for omega_ntr
  real<lower = 0> scale_omega_mtt; // Prior scale parameter for omega_mtt
  
  real<lower = 0> lkj_df_omega;   // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma_p;  // Prior Scale parameter for proportional error
  
  int n_dose;
  array[n_dose] real dosetime;
  array[n_dose] real doseamt;
  
  array[n_subjects] int subj_start_dose;
  array[n_subjects] int subj_end_dose;
  
  int<lower = 0, upper = 1> prior_only;
  
  int<lower = 1, upper = 4> solver; // 1 = rk45, 2 = bdf, 3 = adams, 4 = ckrk
  
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
  
  int n_random = 7;           // Number of random effects
  int n_cmt = 3;              // Absorption, central
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc,
                                      scale_omega_q, scale_omega_vp,
                                      scale_omega_ka, scale_omega_ntr,
                                      scale_omega_mtt}; 
  
  array[n_subjects] int seq_subj = sequence(1, n_subjects); // reduce_sum over subjects
  
  array[n_subjects] real bioav = rep_array(1.0, n_subjects); // Hardcoding, but could be data or a parameter in another situation

  vector[n_cmt] y0 = rep_vector(0.0, n_cmt);
  
}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVQ;
  real<lower = 0> TVVP;
  real<lower = 0.5*(TVCL/TVVC + TVQ/TVVC + TVQ/TVVP +
    sqrt((TVCL/TVVC + TVQ/TVVC + TVQ/TVVP)^2 - 4*TVCL/TVVC*TVQ/TVVP))> TVKA;
  real<lower = 0> TVNTR; 
  real<lower = 0> TVMTT; 
  
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
  vector[n_subjects] eta_ntr;
  vector[n_subjects] eta_mtt;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;
  vector[n_subjects] KA;
  vector[n_subjects] NTR;
  vector[n_subjects] MTT;
  vector[n_subjects] KTR;

  {
  
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP, 
                                                         TVKA, TVNTR, TVMTT});

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
    
    eta_cl = col(eta, 1);
    eta_vc = col(eta, 2);
    eta_q = col(eta, 3);
    eta_vp = col(eta, 4);
    eta_ka = col(eta, 5);
    eta_ntr = col(eta, 6);
    eta_mtt = col(eta, 7);
    CL = col(theta, 1);
    VC = col(theta, 2);
    Q = col(theta, 3);
    VP = col(theta, 4);
    KA = col(theta, 5);
    NTR = col(theta, 6);
    MTT = col(theta, 7);
    KTR = (NTR + 1) ./ MTT;
  
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
  TVNTR ~ lognormal(log(location_tvntr), scale_tvntr); 
  TVMTT ~ lognormal(log(location_tvmtt), scale_tvmtt);

  omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma_p ~ normal(0, scale_sigma_p);
  
  to_vector(Z) ~ std_normal();
  
  // Likelihood
  if(prior_only == 0){
    target += reduce_sum(partial_sum_lupmf, seq_subj, grainsize,
                         dv_obs, dv_obs_id, i_obs,
                         amt, cmt, evid, time, 
                         rate, ii, addl, ss, 
                         subj_start, subj_end, subj_start_dose, subj_end_dose,
                         CL, VC, Q, VP, KA, NTR, KTR,
                         sigma_p,
                         lloq, bloq,
                         n_random, n_subjects, n_total,
                         dosetime, doseamt,
                         bioav, n_cmt, y0, solver);
  }
}
generated quantities{

  real<lower = 0> TVKTR = (TVNTR + 1)/TVMTT;
  
  real<lower = 0> sigma_sq_p = square(sigma_p);
  
  real<lower = 0> omega_cl = omega[1];
  real<lower = 0> omega_vc = omega[2];
  real<lower = 0> omega_q = omega[3];
  real<lower = 0> omega_vp = omega[4];
  real<lower = 0> omega_ka = omega[5];
  real<lower = 0> omega_ntr = omega[6];
  real<lower = 0> omega_mtt = omega[7];

  real<lower = 0> omega_sq_cl = square(omega_cl);
  real<lower = 0> omega_sq_vc = square(omega_vc);
  real<lower = 0> omega_sq_q = square(omega_q);
  real<lower = 0> omega_sq_vp = square(omega_vp);
  real<lower = 0> omega_sq_ka = square(omega_ka);
  real<lower = 0> omega_sq_ntr = square(omega_ntr);
  real<lower = 0> omega_sq_mtt = square(omega_mtt);

  real cor_cl_vc;
  real cor_cl_q;
  real cor_cl_vp;
  real cor_cl_ka;
  real cor_cl_ntr;
  real cor_cl_mtt;
  real cor_vc_q;
  real cor_vc_vp;
  real cor_vc_ka;
  real cor_vc_ntr;
  real cor_vc_mtt;
  real cor_q_vp;
  real cor_q_ka;
  real cor_q_ntr;
  real cor_q_mtt;
  real cor_vp_ka;
  real cor_vp_ntr;
  real cor_vp_mtt;
  real cor_ka_ntr;
  real cor_ka_mtt;
  real cor_ntr_mtt;
  
  real omega_cl_vc;
  real omega_cl_q;
  real omega_cl_vp;
  real omega_cl_ka;
  real omega_cl_ntr;
  real omega_cl_mtt;
  real omega_vc_q;
  real omega_vc_vp;
  real omega_vc_ka;
  real omega_vc_ntr;
  real omega_vc_mtt;
  real omega_q_vp;
  real omega_q_ka;
  real omega_q_ntr;
  real omega_q_mtt;
  real omega_vp_ka;
  real omega_vp_ntr;
  real omega_vp_mtt;
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

    cor_cl_vc = R[1, 2];
    cor_cl_q = R[1, 3];
    cor_cl_vp = R[1, 4];
    cor_cl_ka = R[1, 5];
    cor_cl_ntr = R[1, 6];
    cor_cl_mtt = R[1, 7];
    cor_vc_q = R[2, 3];
    cor_vc_vp = R[2, 4];
    cor_vc_ka = R[2, 5];
    cor_vc_ntr = R[2, 6];
    cor_vc_mtt = R[2, 7];
    cor_q_vp = R[3, 4];
    cor_q_ka = R[3, 5];
    cor_q_ntr = R[3, 6];
    cor_q_mtt = R[3, 7];
    cor_vp_ka = R[4, 5];
    cor_vp_ntr = R[4, 6];
    cor_vp_mtt = R[4, 7];
    cor_ka_ntr = R[5, 6];
    cor_ka_mtt = R[5, 7];
    cor_ntr_mtt = R[6, 7];
    
    omega_cl_vc = Omega[1, 2];
    omega_cl_q = Omega[1, 3];
    omega_cl_vp = Omega[1, 4];
    omega_cl_ka = Omega[1, 5];
    omega_cl_ntr = Omega[1, 6];
    omega_cl_mtt = Omega[1, 7];
    omega_vc_q = Omega[2, 3];
    omega_vc_vp = Omega[2, 4];
    omega_vc_ka = Omega[2, 5];
    omega_vc_ntr = Omega[2, 6];
    omega_vc_mtt = Omega[2, 7];
    omega_q_vp = Omega[3, 4];
    omega_q_ka = Omega[3, 5];
    omega_q_ntr = Omega[3, 6];
    omega_q_mtt = Omega[3, 7];
    omega_vp_ka = Omega[4, 5];
    omega_vp_ntr = Omega[4, 6];
    omega_vp_mtt = Omega[4, 7];
    omega_ka_ntr = Omega[5, 6];
    omega_ka_mtt = Omega[5, 7];
    omega_ntr_mtt = Omega[6, 7];
    
  }

  if(no_gq_predictions == 0){
    
    vector[n_total] dv_ipred;
    array[n_total] vector[n_cmt] x_ipred;
    vector[n_total] dv_pred;
    array[n_total] vector[n_cmt] x_pred;

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
          ode_rk45(transit_2cmt_ode, y0, t0,
                   time[(subj_start[j] + 1):subj_end[j]],
                   CL[j], VC[j], Q[j], VP[j], KA[j], NTR[j], bioav[j], KTR[j],
                   n_dose_subj, dosetime_subj, doseamt_subj);

        x_pred[(subj_start[j] + 1):subj_end[j],] =
          ode_rk45(transit_2cmt_ode, y0, t0,
                   time[(subj_start[j] + 1):subj_end[j]],
                   TVCL, TVVC, TVQ, TVVP, TVKA, TVNTR, bioav[j], TVKTR,
                   n_dose_subj, dosetime_subj, doseamt_subj);
      }else if(solver == 2){
        x_ipred[(subj_start[j] + 1):subj_end[j],] =
          ode_bdf(transit_2cmt_ode, y0, t0,
                   time[(subj_start[j] + 1):subj_end[j]],
                   CL[j], VC[j], Q[j], VP[j], KA[j], NTR[j], bioav[j], KTR[j],
                   n_dose_subj, dosetime_subj, doseamt_subj);

        x_pred[(subj_start[j] + 1):subj_end[j],] =
          ode_bdf(transit_2cmt_ode, y0, t0,
                  time[(subj_start[j] + 1):subj_end[j]],
                  TVCL, TVVC, TVQ, TVVP, TVKA, TVNTR, bioav[j], TVKTR,
                  n_dose_subj, dosetime_subj, doseamt_subj);
      }else if(solver == 3){
        x_ipred[(subj_start[j] + 1):subj_end[j],] =
          ode_adams(transit_2cmt_ode, y0, t0,
                    time[(subj_start[j] + 1):subj_end[j]],
                    CL[j], VC[j], Q[j], VP[j], KA[j], NTR[j], bioav[j], KTR[j],
                    n_dose_subj, dosetime_subj, doseamt_subj);

        x_pred[(subj_start[j] + 1):subj_end[j],] =
          ode_adams(transit_2cmt_ode, y0, t0,
                    time[(subj_start[j] + 1):subj_end[j]],
                    TVCL, TVVC, TVQ, TVVP, TVKA, TVNTR, bioav[j], TVKTR,
                    n_dose_subj, dosetime_subj, doseamt_subj);
      }else{
        x_ipred[(subj_start[j] + 1):subj_end[j],] =
          ode_ckrk(transit_2cmt_ode, y0, t0,
                   time[(subj_start[j] + 1):subj_end[j]],
                   CL[j], VC[j], Q[j], VP[j], KA[j], NTR[j], bioav[j], KTR[j],
                   n_dose_subj, dosetime_subj, doseamt_subj);

        x_pred[(subj_start[j] + 1):subj_end[j],] =
          ode_ckrk(transit_2cmt_ode, y0, t0,
                   time[(subj_start[j] + 1):subj_end[j]],
                   TVCL, TVVC, TVQ, TVVP, TVKA, TVNTR, bioav[j], TVKTR,
                   n_dose_subj, dosetime_subj, doseamt_subj);
      }

      for(k in subj_start[j]:subj_end[j]){
        dv_ipred[k] = fmax(1e-14, x_ipred[k, 2] / VC[j]);
        dv_pred[k] = fmax(1e-14, x_pred[k, 2] / TVVC);
      }

    }

    pred = dv_pred[i_obs];
    ipred = dv_ipred[i_obs];

    res = dv_obs - pred;
    ires = dv_obs - ipred;

    for(i in 1:n_obs){
      real ipred_tmp = ipred[i];
      real sigma_tmp = ipred_tmp*sigma_p;
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
}




