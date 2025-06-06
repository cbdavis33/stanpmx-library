// IV infusion
// Two-compartment PK Model
// IIV on CL, VC, Q, and VP (full covariance matrix)
// exponential error - DV = IPRED*exp(eps)
// Matrix-exponential solution using Torsten 
// Implements threading for within-chain parallelization 
// Deals with BLOQ values by the M3 method (M3 and M4 are equivalent with this
//   error model)
// Covariates: 
//   1) Body Weight on CL, VC, Q, and VP - (wt/70)^theta
//   2) Sex on VC (0 = male, 1 = female) - exp(theta*sex)
//   3) eGFR on CL (continuous) - (eGFR/90)^theta
//   4) Race on VC - exp(theta_vc_{race}*I(race == {race}))

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
  
  real partial_sum_lpmf(array[] int seq_subj, int start, int end,
                        vector dv_obs, array[] int dv_obs_id, array[] int i_obs,
                        array[] real amt, array[] int cmt, array[] int evid, 
                        array[] real time, array[] real rate, array[] real ii, 
                        array[] int addl, array[] int ss,
                        array[] int subj_start, array[] int subj_end, 
                        vector CL, vector VC, vector Q, vector VP,
                        real sigma, 
                        vector lloq, array[] int bloq,
                        int n_random, int n_subjects, int n_total,
                        array[] real bioav, array[] real tlag, int n_cmt){
                           
    real ptarget = 0;
                              
    int N = end - start + 1;    // number of subjects in this slice  
    vector[n_total] dv_ipred;   
    matrix[n_total, n_cmt] x_ipred;
  
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
        
      real ke = CL[j]/VC[j];
      real k_cp = Q[j]/VC[j];
      real k_pc = Q[j]/VP[j];
        
      matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
      K[1, 1] = -(ke + k_cp);
      K[1, 2] = k_pc;
      K[2, 1] = k_cp;
      K[2, 2] = -k_pc;
      
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
        x_ipred[subj_start[j]:subj_end[j], 1] ./ VC[j];
    
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
  
  array[n_subjects] real<lower = 0> wt;                    // baseline bodyweight (kg)
  array[n_subjects] int<lower = 0, upper = 1> sex;         // sex
  array[n_subjects] real<lower = 0> egfr;                  // eGFR
  int<lower = 2> n_races;                                  // number of unique races
  array[n_subjects] int<lower = 1, upper = n_races> race;  // race
  
  real<lower = 0> location_tvcl;  // Prior Location parameter for CL
  real<lower = 0> location_tvvc;  // Prior Location parameter for VC
  real<lower = 0> location_tvq;   // Prior Location parameter for Q
  real<lower = 0> location_tvvp;  // Prior Location parameter for VP
  
  real<lower = 0> scale_tvcl;     // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;     // Prior Scale parameter for VC
  real<lower = 0> scale_tvq;      // Prior Scale parameter for Q
  real<lower = 0> scale_tvvp;     // Prior Scale parameter for VP
  
  real<lower = 0> scale_omega_cl; // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc; // Prior scale parameter for omega_vc
  real<lower = 0> scale_omega_q;  // Prior scale parameter for omega_q
  real<lower = 0> scale_omega_vp; // Prior scale parameter for omega_vp
  
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
  
  int n_random = 4;                    // Number of random effects
  int n_cmt = 2;                       // Number of states in the ODEs
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc, 
                                      scale_omega_q, scale_omega_vp};
  
  array[n_subjects] int seq_subj = linspaced_int_array(n_subjects, 1, n_subjects); // reduce_sum over subjects
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  vector[n_subjects] wt_over_70 = to_vector(wt)/70;
  vector[n_subjects] egfr_over_90 = to_vector(egfr)/90;
  
}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVQ;       
  real<lower = 0> TVVP;
  
  real theta_cl_wt;
  real theta_vc_wt;
  real theta_q_wt;
  real theta_vp_wt;
  real theta_vc_sex;
  real theta_cl_egfr;
  real theta_vc_race2;
  real theta_vc_race3;
  real theta_vc_race4;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
transformed parameters{
  
  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_vc;
  vector[n_subjects] eta_q;
  vector[n_subjects] eta_vp;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;

  {
  
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP});

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    vector[n_races] theta_vc_race = [0, theta_vc_race2, theta_vc_race3, 
                                     theta_vc_race4]';
    
    eta_cl = col(eta, 1);
    eta_vc = col(eta, 2);
    eta_q = col(eta, 3);
    eta_vp = col(eta, 4);
    
    for(j in 1:n_subjects){
      
      real wt_adjustment_cl = wt_over_70[j]^theta_cl_wt;
      real wt_adjustment_vc = wt_over_70[j]^theta_vc_wt;
      real wt_adjustment_q = wt_over_70[j]^theta_q_wt;
      real wt_adjustment_vp = wt_over_70[j]^theta_vp_wt;
      real sex_adjustment_vc = exp(theta_vc_sex*sex[j]);
      real egfr_adjustment_cl = (egfr_over_90[j])^theta_cl_egfr;
      real race_adjustment_vc = exp(theta_vc_race[race[j]]);
      
      row_vector[n_random] theta_j = theta[j]; // access the parameters for subject j
      CL[j] = theta_j[1] * wt_adjustment_cl * egfr_adjustment_cl;
      VC[j] = theta_j[2] * wt_adjustment_vc * race_adjustment_vc * sex_adjustment_vc;
      Q[j] = theta_j[3] * wt_adjustment_q;
      VP[j] = theta_j[4] * wt_adjustment_vp;
      
    }
  
  }
  
}
model{ 
  
  // Priors
  TVCL ~ lognormal(log(location_tvcl), scale_tvcl);
  TVVC ~ lognormal(log(location_tvvc), scale_tvvc);
  TVQ ~ lognormal(log(location_tvq), scale_tvq);
  TVVP ~ lognormal(log(location_tvvp), scale_tvvp);
  
  theta_cl_wt ~ std_normal();
  theta_vc_wt ~ std_normal();
  theta_q_wt ~ std_normal();
  theta_vp_wt ~ std_normal();
  theta_vc_sex ~ std_normal();
  theta_cl_egfr ~ std_normal();
  theta_vc_race2 ~ std_normal();
  theta_vc_race3 ~ std_normal();
  theta_vc_race4 ~ std_normal();

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
                         CL, VC, Q, VP,
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
  real<lower = 0> omega_q = omega[3];
  real<lower = 0> omega_vp = omega[4];

  real<lower = 0> omega_sq_cl = square(omega_cl);
  real<lower = 0> omega_sq_vc = square(omega_vc);
  real<lower = 0> omega_sq_q = square(omega_q);
  real<lower = 0> omega_sq_vp = square(omega_vp);

  real cor_cl_vc;
  real cor_cl_q;
  real cor_cl_vp;
  real cor_vc_q;
  real cor_vc_vp;
  real cor_q_vp;
  real omega_cl_vc;
  real omega_cl_q;
  real omega_cl_vp;
  real omega_vc_q;
  real omega_vc_vp;
  real omega_q_vp;

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
    cor_cl_q = R[1, 3];
    cor_cl_vp = R[1, 4];
    cor_vc_q = R[2, 3];
    cor_vc_vp = R[2, 4];
    cor_q_vp = R[3, 4];
    
    omega_cl_vc = Omega[1, 2];
    omega_cl_q = Omega[1, 3];
    omega_cl_vp = Omega[1, 4];
    omega_vc_q = Omega[2, 3];
    omega_vc_vp = Omega[2, 4];
    omega_q_vp = Omega[3, 4];
    
  }

  if(no_gq_predictions == 0){
    
    vector[n_subjects] CL_new;
    vector[n_subjects] VC_new;
    vector[n_subjects] Q_new;
    vector[n_subjects] VP_new;
    
    vector[n_total] dv_pred;
    matrix[n_total, n_cmt] x_pred;
    vector[n_total] dv_epred;
    matrix[n_total, n_cmt] x_epred;
    vector[n_total] dv_ipred;
    matrix[n_total, n_cmt] x_ipred;
    
    matrix[n_subjects, n_random] eta_new;
    matrix[n_subjects, n_random] theta_new;
    
    vector[n_races] theta_vc_race = [0, theta_vc_race2, theta_vc_race3, 
                                     theta_vc_race4]';
    
    for(i in 1:n_subjects){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(to_row_vector({TVCL, TVVC, TVQ, TVVP}), n_subjects) .* exp(eta_new));

    for(j in 1:n_subjects){
      
      real wt_adjustment_cl = wt_over_70[j]^theta_cl_wt;
      real wt_adjustment_vc = wt_over_70[j]^theta_vc_wt;
      real wt_adjustment_q = wt_over_70[j]^theta_q_wt;
      real wt_adjustment_vp = wt_over_70[j]^theta_vp_wt;
      real sex_adjustment_vc = exp(theta_vc_sex*sex[j]);
      // real egfr_adjustment_cl = (egfr[j]/90)^theta_cl_egfr;
      real egfr_adjustment_cl = (egfr_over_90[j])^theta_cl_egfr;
      real race_adjustment_vc = exp(theta_vc_race[race[j]]);
        
      row_vector[n_random] theta_j_new = theta_new[j]; // access the parameters for subject j's epred
      CL_new[j] = theta_j_new[1] * wt_adjustment_cl * egfr_adjustment_cl;
      VC_new[j] = theta_j_new[2] * wt_adjustment_vc * race_adjustment_vc * sex_adjustment_vc;
      Q_new[j] = theta_j_new[3] * wt_adjustment_q;
      VP_new[j] = theta_j_new[4] * wt_adjustment_vp;
      
      real cl_p = TVCL * wt_adjustment_cl * egfr_adjustment_cl;
      real vc_p = TVVC * wt_adjustment_vc * race_adjustment_vc * sex_adjustment_vc;
      real q_p = TVQ * wt_adjustment_q;
      real vp_p = TVVP * wt_adjustment_vp;
    
      matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
      matrix[n_cmt, n_cmt] K_epred = rep_matrix(0, n_cmt, n_cmt);
      matrix[n_cmt, n_cmt] K_p = rep_matrix(0, n_cmt, n_cmt);
      
      real ke = CL[j]/VC[j];
      real k_cp = Q[j]/VC[j];
      real k_pc = Q[j]/VP[j];
    
      K[1, 1] = -(ke + k_cp);
      K[1, 2] = k_pc;
      K[2, 1] = k_cp;
      K[2, 2] = -k_pc;
      
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
                           
      real ke_new = CL_new[j]/VC_new[j];
      real k_cp_new = Q_new[j]/VC_new[j];
      real k_pc_new = Q_new[j]/VP_new[j];
    
      K_epred[1, 1] = -(ke_new + k_cp_new);
      K_epred[1, 2] = k_pc_new;
      K_epred[2, 1] = k_cp_new;
      K_epred[2, 2] = -k_pc_new;
      
      x_epred[subj_start[j]:subj_end[j], ] =
        pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         K_epred, bioav, tlag)';
      
      real ke_p = cl_p/vc_p;
      real k_cp_p = q_p/vc_p;
      real k_pc_p = q_p/vp_p;
    
      K_p[1, 1] = -(ke_p + k_cp_p);
      K_p[1, 2] = k_pc_p;
      K_p[2, 1] = k_cp_p;
      K_p[2, 2] = -k_pc_p;

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         K_p, bioav, tlag)';
      
      dv_ipred[subj_start[j]:subj_end[j]] =
        x_ipred[subj_start[j]:subj_end[j], 1] ./ VC[j];
        
      dv_epred[subj_start[j]:subj_end[j]] =
        x_epred[subj_start[j]:subj_end[j], 1] ./ VC_new[j];
      
      dv_pred[subj_start[j]:subj_end[j]] =
        x_pred[subj_start[j]:subj_end[j], 1] ./ vc_p;
      
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
