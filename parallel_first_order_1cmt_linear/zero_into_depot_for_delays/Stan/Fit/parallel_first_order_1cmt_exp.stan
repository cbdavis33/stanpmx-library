// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// This is flexible so that any number of parallel absorption processes can take 
//   place (n_depots >= 2).
// A delay in absorption for each process is implemented with a zero-order 
//   distributive delay (like an infusion into the depot) to bring about a
//   delayed absorption. The user can choose whether the fastest process has a 
//   delay by setting n_depots_with_delay to (n_depots - 1) for no delay or to
//   n_depots for a delay in the data input
// There will be n_depots KA values for each subject. i.e., TVKA and KAi are 
//   vectors. KA is intended to be ordered so that 
//   KA_1 < KA_2 < ... < KA_n_depots. This doesn't do this exactly, but it will 
//   have TVKA_1 < TVKA_2 < ... < TVKA_n_depots and hope the individual effects 
//   don't mess that up
// There will be n_depots_with_delay DUR values for each subject. i.e., TVDUR 
//   and DURi are vectors. DUR is intended to be ordered so that 
//   DUR_1 > DUR_2 > ... > DUR_n_depots_with_delay. If 
//   n_depots_with_delay = n_depots, then the fastest process has no delay. This
//   doesn't do this exactly, but it will have 
//   TVDUR_1 > TVDUR_2 > ... > TVDUR_(n_depots - 1) and hope the individual
//   effects don't mess that up
// The prior on TVKA is a half-normal(0, scale_tvka). scale_tvka can be an array
//   of length 1 (which gives the same prior to each TVKA_i) OR it can be an 
//   array of length n_depots. This is meant to provide maximum flexibility. The
//   user must provide length_tvka_prior. This MUST be either 1 or n_depots
// The prior on TVDUR is a half-normal(0, scale_tvka). scale_tvdur can be an 
//   array of length 1 (which gives the same prior to each TVDUR_i) OR it can be  
//   an array of length n_depots_with_delay. This is meant to provide maximum 
//   flexibility. The user must provide length_tvdur_prior. This MUST be either 
//   1 or n_depots_with_delay
// TVFRAC is a simplex of length n_depots that tells how much of the dose goes 
//   into each absorption process. There is no IIV on FRAC
// IIV on CL, VC, KA, DUR (full covariance matrix)
// Matrix-exponential solution using Torsten
// exponential error - DV = IPRED*exp(eps)
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
                        array[] real time, array[] real ii,  
                        array[] int addl, array[] int ss,
                        array[] int subj_start, array[] int subj_end, 
                        vector CL, vector VC, 
                        array[] row_vector KA, array[] row_vector DUR,
                        real sigma, 
                        vector lloq, array[] int bloq,
                        int n_random, int n_subjects, int n_total,
                        array[] real bioav, array[] real tlag, int n_cmt, 
                        int n_depots, int n_depots_with_delay,
                        data array[,] real x_r, data array[,] int x_i){
                           
    real ptarget = 0;
                              
    int N = end - start + 1;    // number of subjects in this slice  
    vector[n_total] dv_ipred;   
    matrix[n_total, n_cmt] x_ipred;
    
    array[n_total] real rate;
  
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
      
      for(i in subj_start[j]:subj_end[j]){

        if(cmt[i] <= n_depots_with_delay){
          rate[i] = amt[i]/DUR[j, cmt[i]];
        }else{
          rate[i] = 0;
        }
        if(is_inf(rate[i])) rate[i] = 0;
      
      }

      matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
        
      for(i in 1:n_depots){
        K[i, i] = -KA[j, i];
        K[(n_depots + 1), i] = KA[j, i];
      }
        
      K[(n_depots + 1), (n_depots + 1)] = -CL[j]/VC[j];
        
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
                           
      dv_ipred[subj_start[j]:subj_end[j]] = 
                    x_ipred[subj_start[j]:subj_end[j], (n_depots + 1)] ./ VC[j];
    
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
  array[n_total] real ii;
  array[n_total] int addl;
  array[n_total] int ss;
  array[n_total] real time;
  vector<lower = 0>[n_total] dv;
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;
  vector[n_total] lloq;
  array[n_total] int bloq;
  
  int<lower = 2> n_depots; 
  int<lower = n_depots - 1, upper = n_depots> n_depots_with_delay;
  
  int<lower = 1, upper = n_depots> length_tvka_prior; 
  int<lower = 1, upper = n_depots_with_delay> length_tvdur_prior;
  
  real<lower = 0> location_tvcl;  // Prior Location parameter for CL
  real<lower = 0> location_tvvc;  // Prior Location parameter for VC
  
  array[length_tvka_prior] real<lower = 0> scale_tvka;   // Prior Scale parameter for KA
  real<lower = 0> scale_tvcl;                            // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;                            // Prior Scale parameter for VC
  array[length_tvdur_prior] real<lower = 0> scale_tvdur; // Prior Scale parameter for DUR
  
  array[n_depots] real<lower = 0> scale_omega_ka;             // Prior scale parameter for omega_ka
  real<lower = 0> scale_omega_cl;                             // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc;                             // Prior scale parameter for omega_vc
  array[n_depots_with_delay] real<lower = 0> scale_omega_dur; // Prior scale parameter for omega_dur
  
  real<lower = 0> alpha_tvfrac; // Prior for dirichlet distribution for TVFRAC
  
  real<lower = 0> lkj_df_omega;   // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma;    // Prior Scale parameter for exponential error
  
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
  
  int n_random = n_depots + n_depots_with_delay + 2; // Number of random effects
  int n_cmt = n_depots + 1;                          // Number of states in the ODEs
  
  array[n_random] real scale_omega = 
    append_array(append_array(scale_omega_ka, {scale_omega_cl, scale_omega_vc}), 
                 scale_omega_dur); 
  
  array[n_subjects] int seq_subj = sequence(1, n_subjects); // reduce_sum over subjects
  
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  array[1, 2] real x_r = {{8675309, 5555555}}; // This is a placeholder of nonsense so I can put it in the Torsten function
  
  array[1, 1] int x_i = {{n_depots}};
  
}
parameters{ 
  
  positive_ordered[n_depots] TVKA;
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  positive_ordered[n_depots_with_delay] TVDUR_backward;
  
  simplex[n_depots] TVFRAC;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
transformed parameters{
  
  array[n_subjects] row_vector[n_depots] eta_ka;
  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_vc;
  array[n_subjects] row_vector[n_depots_with_delay] eta_dur;
  
  array[n_subjects] row_vector[n_depots] KA;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  array[n_subjects] row_vector[n_depots_with_delay] DUR;
  vector[n_subjects] KE;
  
  vector[n_depots_with_delay] TVDUR = reverse(TVDUR_backward);

  {
    
    row_vector[n_random] typical_values = 
                            append_col(append_col(TVKA', [TVCL, TVVC]), TVDUR');

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
    
    eta_cl = col(eta, n_depots + 1);
    eta_vc = col(eta, n_depots + 2);
    
    CL = col(theta, n_depots + 1);
    VC = col(theta, n_depots + 2);

    KE = CL ./ VC;
    
    for(j in 1:n_subjects){
      
      eta_ka[j] = eta[j, 1:n_depots];
      eta_dur[j] = eta[j, (n_depots + 2 + 1):n_random];
      
      KA[j] = theta[j, 1:n_depots];
      DUR[j] = theta[j, (n_depots + 2 + 1):n_random];
      
    }
  
  }
  
}
model{ 
  
  // Priors
  if(length_tvka_prior == 1){
    TVKA ~ normal(0, scale_tvka[1]);
  }else{
    TVKA ~ normal(0, scale_tvka);
  }
  TVCL ~ lognormal(log(location_tvcl), scale_tvcl);
  TVVC ~ lognormal(log(location_tvvc), scale_tvvc);
  if(length_tvdur_prior == 1){
    TVDUR_backward ~ normal(0, scale_tvdur[1]);
  }else{
    TVDUR_backward ~ normal(0, scale_tvdur);
  }
  
  TVFRAC ~ dirichlet(rep_vector(alpha_tvfrac, n_depots));

  omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma ~ normal(0, scale_sigma);
  
  to_vector(Z) ~ std_normal();
  
  // Likelihood
  if(prior_only == 0){
    target += reduce_sum(partial_sum_lupmf, seq_subj, grainsize,
                         dv_obs, dv_obs_id, i_obs,
                         amt, cmt, evid, time,
                         ii, addl, ss, subj_start, subj_end, 
                         CL, VC, KA, DUR,
                         sigma,
                         lloq, bloq,
                         n_random, n_subjects, n_total,
                         append_array(to_array_1d(TVFRAC), {1.0}), 
                         tlag, n_cmt, n_depots, n_depots_with_delay,
                         x_r, x_i);
  }
}
generated quantities{
  
  real<lower = 0> sigma_sq = square(sigma);
  
  vector<lower = 0>[n_depots] omega_ka = omega[1:n_depots];
  real<lower = 0> omega_cl = omega[n_depots + 1];
  real<lower = 0> omega_vc = omega[n_depots + 2];
  vector<lower = 0>[n_depots_with_delay] omega_dur = omega[(n_depots + 2 + 1):n_random];
  
  vector<lower = 0>[n_depots] omega_sq_ka = square(omega_ka);
  real<lower = 0> omega_sq_cl = square(omega_cl);
  real<lower = 0> omega_sq_vc = square(omega_vc);
  vector<lower = 0>[n_depots_with_delay] omega_sq_dur = square(omega_dur);

  // Normally this is where I'd try to return named correlations and covariances
  // (e.g. cor_cl_ka, omega_cl_ka), but since n_depots can vary, there's no 
  // good way to return them. It'll just have to be the whole correlation matrix
  // and whole covariance matrix, and the user should sort it out in R
  
  matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
  matrix[n_random, n_random] Omega = quad_form_diag(R, omega);
  
  vector[n_obs] ipred;
  vector[n_obs] pred;
  vector[n_obs] dv_ppc;
  vector[n_obs] log_lik;
  vector[n_obs] res;
  vector[n_obs] wres;
  vector[n_obs] ires;
  vector[n_obs] iwres;
 
  if(no_gq_predictions == 0){

    vector[n_total] dv_pred;
    matrix[n_total, n_cmt] x_pred;
    vector[n_total] dv_ipred;
    matrix[n_total, n_cmt] x_ipred;
    
    array[n_total] real rate;
    array[n_total] real rate_p;
    
    array[n_cmt] real bioav = append_array(to_array_1d(TVFRAC), {1.0}); 

    for(j in 1:n_subjects){
      
      // The x_p values are here to make it easier to incorporate covariates later
      vector[n_depots] ka_p = TVKA;
      real cl_p = TVCL;
      real vc_p = TVVC;
      vector[n_depots_with_delay] dur_p = TVDUR;
                    
      for(i in subj_start[j]:subj_end[j]){

        if(cmt[i] <= n_depots_with_delay){
          
          rate[i] = amt[i]/DUR[j, cmt[i]];
          rate_p[i] = amt[i]/dur_p[cmt[i]];
          
        }else{
          rate[i] = 0;
          rate_p[i] = 0;
        }
        if(is_inf(rate[i])) rate[i] = 0;
        if(is_inf(rate_p[i])) rate_p[i] = 0;
      
      }
    
      matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
      matrix[n_cmt, n_cmt] K_tv = rep_matrix(0, n_cmt, n_cmt);
        
      for(i in 1:n_depots){
        K[i, i] = -KA[j, i];
        K[(n_depots + 1), i] = KA[j, i];
      }
        
      K[(n_depots + 1), (n_depots + 1)] = -KE[j];
      
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
      
        
      for(i in 1:n_depots){
        K_tv[i, i] = -ka_p[i];
        K_tv[(n_depots + 1), i] = ka_p[i];
      }
      
      K_tv[(n_depots + 1), (n_depots + 1)] = -cl_p/vc_p;
        x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate_p[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         K_tv, bioav, tlag)';
                         
      dv_ipred[subj_start[j]:subj_end[j]] =
        x_ipred[subj_start[j]:subj_end[j], (n_depots + 1)] ./ VC[j];
      
      dv_pred[subj_start[j]:subj_end[j]] =
        x_pred[subj_start[j]:subj_end[j], (n_depots + 1)] ./ vc_p;
      
    }

    pred = dv_pred[i_obs];
    ipred = dv_ipred[i_obs];

    res = log(dv_obs) - log(pred);
    ires = log(dv_obs) - log(ipred);

    for(i in 1:n_obs){
      real log_ipred_tmp = log(ipred[i]);
      dv_ppc[i] = lognormal_rng(log_ipred_tmp, sigma);
      if(bloq_obs[i] == 1){
        log_lik[i] = lognormal_lcdf(lloq_obs[i] | log_ipred_tmp, sigma);
      }else{
        log_lik[i] = lognormal_lpdf(dv_obs[i] | log_ipred_tmp, sigma);
      }
      wres[i] = res[i]/sigma;
      iwres[i] = ires[i]/sigma;
    }
  }
}

