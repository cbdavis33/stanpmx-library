// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// Nonlinear Bioavailability - Pin the bioavailability to be 1 at
//   'reference_dose' mg. Bioavailability decreases in the form 
//   bioav = 1 - FMAX*(dose - reference_dose)^HILL/(F50^HILL + (dose - reference_dose)^HILL)
// IIV on CL, VC, KA (full covariance matrix)
// exponential error - DV = IPRED*exp(eps)
// Matrix-exponential solution using Torsten (the matrix-exponential seems to be
//   faster than the analytical solution for this model)
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
  
  vector depot_1cmt_ode(real t, vector y, array[] real params, 
                        array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real ke = cl/vc;
    
    vector[2] dydt;

    dydt[1] = -ka*y[1];           // depot
    dydt[2] = ka*y[1] - ke*y[2];  // central
    
    return dydt;
  }
  
  real partial_sum_lpmf(array[] int seq_subj, int start, int end,
                        vector dv_obs, array[] int dv_obs_id, array[] int i_obs,
                        array[] real amt, array[] int cmt, array[] int evid, 
                        array[] real time, array[] real rate, array[] real ii, 
                        array[] int addl, array[] int ss,
                        array[] int subj_start, array[] int subj_end, 
                        vector CL, vector VC, vector KA, vector BIOAV,
                        real sigma, 
                        vector lloq, array[] int bloq,
                        int n_random, int n_subjects, int n_total,
                        array[] real tlag, int n_cmt){
                           
    real ptarget = 0;
                              
    int N = end - start + 1;    // number of subjects in this slice  
    vector[n_total] dv_ipred;   
    matrix[n_total, 2] x_ipred;
  
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
      
      matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
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
                         K, {BIOAV[j], 1}, tlag)';
                      
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
  real<lower = 0> location_tvka;    // Prior Location parameter for KA
  real<lower = 0, upper = 1> location_fmax; // Prior location parameter for FMAX (0.5)
  real<lower = 0> location_f50;    // Prior location parameter for F50 (amount above reference_dose)
  real<lower = 0> location_hill;   // Prior location parameter for hill parameter
  
  real<lower = 0> scale_tvcl;     // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;     // Prior Scale parameter for VC
  real<lower = 0> scale_tvka;     // Prior Scale parameter for KA
  real<lower = 0> scale_fmax;     // Prior Scale parameter for FMAX (3.5, kappa in beta_proportion())
  real<lower = 0> scale_f50;      // Prior Scale parameter for F50
  real<lower = 0> scale_hill;     // Prior Scale parameter for hill parameter
  
  real<lower = 0> scale_omega_cl;    // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc;    // Prior scale parameter for omega_vc
  real<lower = 0> scale_omega_ka;    // Prior scale parameter for omega_ka
  
  real<lower = 0> reference_dose;  // reference dose for bioavailability (lowest dose)
  array[n_subjects] real<lower = reference_dose> dose; // dose in mg for each subject
  
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
  
  int n_random = 3;                    // Number of random effects
  int n_cmt = 2;                       // Number of states in the ODEs
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc, 
                                      scale_omega_ka}; 
  
  array[n_subjects] int seq_subj = sequence(1, n_subjects); // reduce_sum over subjects
  
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  vector[n_subjects] dose_minus_reference = to_vector(dose) - reference_dose;
  
}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = TVCL/TVVC> TVKA;
  
  real<lower = 0, upper = 1> FMAX;
  real<lower = 0> F50;
  real<lower = 0> HILL;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
transformed parameters{
  
  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_vc;
  vector[n_subjects] eta_ka;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  vector[n_subjects] BIOAV;
  vector[n_subjects] KE;

  {
  
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA});

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
    
    eta_cl = col(eta, 1);
    eta_vc = col(eta, 2);
    eta_ka = col(eta, 3);
    CL = col(theta, 1);
    VC = col(theta, 2);
    KA = col(theta, 3);
    BIOAV = 1 - (FMAX*dose_minus_reference^HILL./(F50^HILL + dose_minus_reference^HILL));
    KE = CL ./ VC;
  
  }
  
}
model{ 
  
  // Priors
  TVCL ~ lognormal(log(location_tvcl), scale_tvcl);
  TVVC ~ lognormal(log(location_tvvc), scale_tvvc);
  TVKA ~ lognormal(log(location_tvka), scale_tvka) T[TVCL/TVVC, ];
  
  FMAX ~ beta_proportion(location_fmax, scale_fmax); 
  F50 ~ lognormal(log(location_f50), scale_f50);
  HILL ~ lognormal(log(location_hill), scale_hill);
  
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
                         CL, VC, KA, BIOAV,
                         sigma,
                         lloq, bloq,
                         n_random, n_subjects, n_total,
                         tlag, n_cmt);
  }
}
generated quantities{
  
  real<lower = 0> sigma_sq = square(sigma);

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
    cor_cl_ka = R[1, 3];
    cor_vc_ka = R[2, 3];
    
    omega_cl_vc = Omega[1, 2];
    omega_cl_ka = Omega[1, 3];
    omega_vc_ka = Omega[2, 3];
    
  }

  if(no_gq_predictions == 0){
    
    vector[n_total] dv_pred;
    matrix[n_total, 2] x_pred;
    vector[n_total] dv_ipred;
    matrix[n_total, 2] x_ipred;
    
    for(j in 1:n_subjects){
      
      real bioav_p = BIOAV[j];
      
      matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
      matrix[n_cmt, n_cmt] K_tv = rep_matrix(0, n_cmt, n_cmt);
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
                         K, {BIOAV[j], 1}, tlag)';
                           
      K_tv[1, 1] = -TVKA;
      K_tv[2, 1] = TVKA;
      K_tv[2, 2] = -TVCL/TVVC;

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         K_tv, {bioav_p, 1}, tlag)';
      
      dv_ipred[subj_start[j]:subj_end[j]] =
        x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
      
      dv_pred[subj_start[j]:subj_end[j]] =
        x_pred[subj_start[j]:subj_end[j], 2] ./ TVVC;
      
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


