// IV infusion
// One-compartment PK Model
// IIV on CL and VC (full covariance matrix)
// proportional error - DV = IPRED*(1 + eps_p)
// Any of analytical, matrix-exponential, or general ODE solution (rk45 or bdf) 
//   using Torsten (user's choice)
// Implements threading for within-chain parallelization 
// Deals with BLOQ values by the "CDF trick" (M4)
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   truncates the likelihood below at 0
// For PPC, it generates values from a normal that is truncated below at 0

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
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y; 

  }
  
  vector iv_1cmt_ode(real t, vector y, array[] real params, 
                     array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    
    real ke = cl/vc;
    
    vector[1] dydt;

    dydt[1] = -ke*y[1];  // central
    
    return dydt;
  }
  
  real partial_sum_lpmf(array[] int seq_subj, int start, int end,
                        vector dv_obs, array[] int dv_obs_id, array[] int i_obs,
                        array[] real amt, array[] int cmt, array[] int evid, 
                        array[] real time, array[] real rate, array[] real ii, 
                        array[] int addl, array[] int ss,
                        array[] int subj_start, array[] int subj_end, 
                        vector CL, vector VC, 
                        real sigma_p, 
                        vector lloq, array[] int bloq,
                        int n_random, int n_subjects, int n_total,
                        array[] real bioav, array[] real tlag, int n_cmt, 
                        int solver){
                           
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
      
      if(solver == 1){
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           {CL[j], VC[j], 0})';
          
      }else if(solver == 2){
        
        real ke = CL[j]/VC[j];
        
        matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
        K[1, 1] = -ke;
      
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
                           
      }else if(solver == 3){
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(iv_1cmt_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j]}, bioav, tlag)';
                           
      }else{
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_bdf(iv_1cmt_ode,
                        n_cmt,
                        time[subj_start[j]:subj_end[j]],
                        amt[subj_start[j]:subj_end[j]],
                        rate[subj_start[j]:subj_end[j]],
                        ii[subj_start[j]:subj_end[j]],
                        evid[subj_start[j]:subj_end[j]],
                        cmt[subj_start[j]:subj_end[j]],
                        addl[subj_start[j]:subj_end[j]],
                        ss[subj_start[j]:subj_end[j]],
                        {CL[j], VC[j]}, bioav, tlag)';
                         
      }
                      
      dv_ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], n_cmt] ./ VC[j];
    
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
  
  real<lower = 0> location_tvcl;  // Prior Location parameter for CL
  real<lower = 0> location_tvvc;  // Prior Location parameter for VC
  
  real<lower = 0> scale_tvcl;     // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;     // Prior Scale parameter for VC
  
  real<lower = 0> scale_omega_cl; // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc; // Prior scale parameter for omega_vc
  
  real<lower = 0> lkj_df_omega;   // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma_p;  // Prior Scale parameter for proportional error
  
  int<lower = 0, upper = 1> prior_only; // Want to simulate from the prior?
  int<lower = 0, upper = prior_only> no_gq_predictions; // Leave out PREDS and IPREDS in 
                                                        // generated quantities. Useful
                                                        // for simulating prior parameters
                                                        // but don't want prior predictions
  
  int<lower = 1, upper = 4> solver; // 1 = analytical, 2 = matrix exponential, 3 = rk45, 4 = bdf
 
}
transformed data{ 
  
  int grainsize = 1;
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 2;                    // Number of random effects
  int n_cmt = (solver == 1) ? 2 : 1;   // Number of states in the ODEs
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc}; 
  
  array[n_subjects] int seq_subj = linspaced_int_array(n_subjects, 1, n_subjects); // reduce_sum over subjects
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
}
transformed parameters{
  
  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_vc;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KE;

  {
  
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC});

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
    
    eta_cl = col(eta, 1);
    eta_vc = col(eta, 2);
    CL = col(theta, 1);
    VC = col(theta, 2);
    KE = CL ./ VC;
  
  }
  
}
model{ 
  
  // Priors
  TVCL ~ lognormal(log(location_tvcl), scale_tvcl);
  TVVC ~ lognormal(log(location_tvvc), scale_tvvc);

  omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma_p ~ normal(0, scale_sigma_p);
  
  to_vector(Z) ~ std_normal();
  
  // Likelihood
  if(prior_only == 0){
    target += reduce_sum(partial_sum_lupmf, seq_subj, grainsize,
                         dv_obs, dv_obs_id, i_obs,
                         amt, cmt, evid, time, 
                         rate, ii, addl, ss, subj_start, subj_end, 
                         CL, VC,
                         sigma_p,
                         lloq, bloq,
                         n_random, n_subjects, n_total,
                         bioav, tlag, n_cmt, solver);
  }
}
generated quantities{
  
  real<lower = 0> sigma_sq_p = square(sigma_p);

  real<lower = 0> omega_cl = omega[1];
  real<lower = 0> omega_vc = omega[2];

  real<lower = 0> omega_sq_cl = square(omega_cl);
  real<lower = 0> omega_sq_vc = square(omega_vc);

  real cor_cl_vc;
  real omega_cl_vc;

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
    omega_cl_vc = Omega[1, 2];
    
  }
 
  if(no_gq_predictions == 0){

    vector[n_subjects] CL_new;
    vector[n_subjects] VC_new;
    vector[n_subjects] KE_new;
    
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
    theta_new = (rep_matrix(to_row_vector({TVCL, TVVC}), n_subjects) .* exp(eta_new));

    for(j in 1:n_subjects){
      
      row_vector[n_random] theta_j_new = theta_new[j]; // access the parameters for subject j's epred
      
      real cl_p = TVCL;
      real vc_p = TVVC;
      
      CL_new[j] = theta_j_new[1];
      VC_new[j] = theta_j_new[2];
      KE_new[j] = CL_new[j]/VC_new[j];
      
      if(solver == 1){
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           {CL[j], VC[j], 0})';
                           
        x_epred[subj_start[j]:subj_end[j],] =
          pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           {CL_new[j], VC_new[j], 0})';
                           
        x_pred[subj_start[j]:subj_end[j],] =
          pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           {cl_p, vc_p, 0})';
          
      }else if(solver == 2){
        
        matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
        matrix[n_cmt, n_cmt] K_epred = rep_matrix(0, n_cmt, n_cmt);
        matrix[n_cmt, n_cmt] K_p = rep_matrix(0, n_cmt, n_cmt);
      
        K[1, 1] = -KE[j];
      
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
                             
        K_epred[1, 1] = -KE_new[j];
      
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
      
        K_p[1, 1] = -cl_p/vc_p;

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
                           
      }else if(solver == 3){
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(iv_1cmt_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j]}, bioav, tlag)';
                         
        x_epred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(iv_1cmt_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL_new[j], VC_new[j], 0}, bioav, tlag)';
                         
        x_pred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(iv_1cmt_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {cl_p, vc_p, 0}, bioav, tlag)';
                           
      }else{
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_bdf(iv_1cmt_ode,
                        n_cmt,
                        time[subj_start[j]:subj_end[j]],
                        amt[subj_start[j]:subj_end[j]],
                        rate[subj_start[j]:subj_end[j]],
                        ii[subj_start[j]:subj_end[j]],
                        evid[subj_start[j]:subj_end[j]],
                        cmt[subj_start[j]:subj_end[j]],
                        addl[subj_start[j]:subj_end[j]],
                        ss[subj_start[j]:subj_end[j]],
                        {CL[j], VC[j]}, bioav, tlag)';
                         
        x_epred[subj_start[j]:subj_end[j],] =
          pmx_solve_bdf(iv_1cmt_ode,
                        n_cmt,
                        time[subj_start[j]:subj_end[j]],
                        amt[subj_start[j]:subj_end[j]],
                        rate[subj_start[j]:subj_end[j]],
                        ii[subj_start[j]:subj_end[j]],
                        evid[subj_start[j]:subj_end[j]],
                        cmt[subj_start[j]:subj_end[j]],
                        addl[subj_start[j]:subj_end[j]],
                        ss[subj_start[j]:subj_end[j]],
                        {CL_new[j], VC_new[j], 0}, bioav, tlag)';
                         
        x_pred[subj_start[j]:subj_end[j],] =
          pmx_solve_bdf(iv_1cmt_ode,
                        n_cmt,
                        time[subj_start[j]:subj_end[j]],
                        amt[subj_start[j]:subj_end[j]],
                        rate[subj_start[j]:subj_end[j]],
                        ii[subj_start[j]:subj_end[j]],
                        evid[subj_start[j]:subj_end[j]],
                        cmt[subj_start[j]:subj_end[j]],
                        addl[subj_start[j]:subj_end[j]],
                        ss[subj_start[j]:subj_end[j]],
                        {cl_p, vc_p, 0}, bioav, tlag)';       
                         
      }
      
      dv_ipred[subj_start[j]:subj_end[j]] =
        x_ipred[subj_start[j]:subj_end[j], n_cmt] ./ VC[j];
        
      dv_epred[subj_start[j]:subj_end[j]] =
        x_epred[subj_start[j]:subj_end[j], n_cmt] ./ VC_new[j];
      
      dv_pred[subj_start[j]:subj_end[j]] =
        x_pred[subj_start[j]:subj_end[j], n_cmt] ./ vc_p;
      
    }

    pred = dv_pred[i_obs];
    epred_stan = dv_epred[i_obs];
    ipred = dv_ipred[i_obs];

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
