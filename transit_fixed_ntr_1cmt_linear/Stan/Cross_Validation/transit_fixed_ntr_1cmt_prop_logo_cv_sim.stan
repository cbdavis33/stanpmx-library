// One-compartment PK Model with Transit Compartment absorption
// IIV on CL, VC, KA, MTT (full covariance matrix)
// n_transit is fixed to a positive integer that can be defined as data
// The 0th transit compartment is where the dosing happens (cmt = 1). This could
//   also be called the Depot
// proportional error - DV = IPRED*(1 + eps_p)
// Matrix-exponential solution using Torsten
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   truncates the likelihood below at 0
// n_sim has to be extremely high (~ 5000-10000 for this model) for PSIS-LOO to
//   be useful

data{
  
  int n_subjects;
  int n_total;
  int n_obs;
  array[n_obs] int i_obs;
  array[n_total] int ID;
  array[n_total] real amt;
  array[n_total] int cmt;
  array[n_total] int evid;
  array[n_total] int rate;
  array[n_total] real ii;
  array[n_total] int addl;
  array[n_total] int ss;
  array[n_total] real time;
  vector<lower = 0>[n_total] dv;
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;
  vector[n_total] lloq;
  array[n_total] int bloq;
  
  int<lower = 1> n_sim;
  
  int<lower = 1> n_transit;
 
}
transformed data{ 
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 4;           // Number of random effects
  int n_cmt = n_transit + 3;  // Depot, tr_1, ..., tr_n, absorption, central
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt);
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  real log_n_sim = log(n_sim); // This only needs calculated once for efficiency
                               // Probably won't make much difference

}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = TVCL/TVVC> TVKA;
  real<lower = 0> TVMTT;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{
 
  vector[n_subjects] log_lik_subj;
 
  {
    
    vector[n_random] typical_values = to_vector({TVCL, TVVC, TVKA, TVMTT});
    array[n_subjects] row_vector[n_sim] log_lik_subj_sim = 
      rep_array(rep_row_vector(0, n_sim), n_subjects);
    
    for(m in 1:n_sim){
      
      matrix[n_total, n_cmt] x_epred;
      vector[n_total] dv_epred;
      vector[n_obs] epred_stan;
      
      vector[n_random] eta_new = 
            multi_normal_cholesky_rng(rep_vector(0, n_random),
                                      diag_pre_multiply(omega, L)); 
      vector[n_random] theta_new = typical_values .* exp(eta_new);
      
      for(j in 1:n_subjects){
        
        // // placeholder for making adjustments for covariates if necessary
        // real wt_over_70 = wt[j]/70;
        // real wt_adjustment_cl = wt_over_70^theta_cl_wt;
        // real wt_adjustment_vc = wt_over_70^theta_vc_wt;
        
        real CL_new = theta_new[1]; // * weight_adjustment_cl;
        real VC_new = theta_new[2]; // * weight_adjustment_vc;
        real KA_new = theta_new[3];
        real MTT_new = theta_new[4];
        real KTR_new = (n_transit + 1)/MTT_new;
                           
        matrix[n_cmt, n_cmt] K_epred = rep_matrix(0, n_cmt, n_cmt);
        
        for(i in 1:(n_transit + 1)){
          K_epred[i, i] = -KTR_new;
          K_epred[(i + 1), i] = KTR_new;
        }
        
        K_epred[(n_transit + 2), (n_transit + 2)] = -KA_new;
        K_epred[(n_transit + 3), (n_transit + 2)] = KA_new;
        K_epred[(n_transit + 3), (n_transit + 3)] = -CL_new/VC_new;
      
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
                           
        dv_epred[subj_start[j]:subj_end[j]] = 
          x_epred[subj_start[j]:subj_end[j], (n_transit + 3)] ./ VC_new;
        
      }
      
      epred_stan = dv_epred[i_obs];
      
      for(i in 1:n_obs){
        real epred_tmp = epred_stan[i];
        real sigma_tmp = epred_tmp*sigma_p;
        if(bloq_obs[i] == 1){
          log_lik_subj_sim[dv_obs_id[i], m] += 
            log_diff_exp(normal_lcdf(lloq_obs[i] | epred_tmp, sigma_tmp),
                         normal_lcdf(0.0 | epred_tmp, sigma_tmp)) -
            normal_lccdf(0.0 | epred_tmp, sigma_tmp);
        }else{
          log_lik_subj_sim[dv_obs_id[i], m] += 
            normal_lpdf(dv_obs[i] | epred_tmp, sigma_tmp) -
                        normal_lccdf(0.0 | epred_tmp, sigma_tmp);
        }

      }
     
    }
    
    for(j in 1:n_subjects){
      log_lik_subj[j] = log_sum_exp(log_lik_subj_sim[j]) - log_n_sim;
    }
    
  }
}

