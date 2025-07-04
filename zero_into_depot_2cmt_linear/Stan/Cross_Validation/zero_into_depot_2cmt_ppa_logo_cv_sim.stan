// First Order Absorption (oral/subcutaneous)
// Two-compartment PK Model
// Zero-order distributive delay (like an infusion into the gut) to bring about 
//   a delayed absorption
// IIV on CL, VC, Q, VP, KA, DUR (full covariance matrix)
// proportional plus additive error - DV = IPRED*(1 + eps_p) + eps_a
// Analytical solution
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
 
}
transformed data{ 
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 6; // Number of random effects
  int n_cmt = 3;    // Number of compartments - depot, central, peripheral
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt);
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  real log_n_sim = log(n_sim); // This only needs calculated once for efficiency
                               // Probably won't make much difference

}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVQ;       
  real<lower = 0> TVVP;
  real<lower = 0.5*(TVCL/TVVC + TVQ/TVVC + TVQ/TVVP +
        sqrt((TVCL/TVVC + TVQ/TVVC + TVQ/TVVP)^2 - 4*TVCL/TVVC*TVQ/TVVP))> TVKA;
  real<lower = 0> TVDUR; 
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  vector<lower = 0>[2] sigma;
  cholesky_factor_corr[2] L_Sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{
 
  vector[n_subjects] log_lik_subj;
 
  {
    
    vector[n_random] typical_values = to_vector({TVCL, TVVC, TVQ, TVVP, TVKA, 
                                                 TVDUR});
    array[n_subjects] row_vector[n_sim] log_lik_subj_sim = 
      rep_array(rep_row_vector(0, n_sim), n_subjects);
      
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
    
    real sigma_sq_p = Sigma[1, 1];
    real sigma_sq_a = Sigma[2, 2];
    real sigma_p_a = Sigma[1, 2];
    
    for(m in 1:n_sim){
      
      matrix[n_total, n_cmt] x_epred;
      vector[n_total] dv_epred;
      vector[n_obs] epred_stan;
      
      vector[n_random] eta_new = 
            multi_normal_cholesky_rng(rep_vector(0, n_random),
                                      diag_pre_multiply(omega, L)); 
      vector[n_random] theta_new = typical_values .* exp(eta_new);
      
      array[n_total] real rate;
      
      for(j in 1:n_subjects){
        
        // // placeholder for making adjustments for covariates if necessary
        // real wt_over_70 = wt[j]/70;
        // real wt_adjustment_cl = wt_over_70^theta_cl_wt;
        // real wt_adjustment_vc = wt_over_70^theta_vc_wt;
        // real wt_adjustment_q = wt_over_70^theta_q_wt;
        // real wt_adjustment_vp = wt_over_70^theta_vp_wt;
        
        real CL_new = theta_new[1]; // * weight_adjustment_cl;
        real VC_new = theta_new[2]; // * weight_adjustment_vc;
        real Q_new = theta_new[3];  // * weight_adjustment_q;
        real VP_new = theta_new[4]; // * weight_adjustment_vp;
        real KA_new = theta_new[5];
        real DUR_new = theta_new[6];
        
        rate[subj_start[j]:subj_end[j]] =  
           to_array_1d(to_vector(amt[subj_start[j]:subj_end[j]]) ./ DUR_new);
        
        x_epred[subj_start[j]:subj_end[j],] =
          pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           {CL_new, Q_new, VC_new, VP_new, KA_new})';
                           
        dv_epred[subj_start[j]:subj_end[j]] = 
          x_epred[subj_start[j]:subj_end[j], 2] ./ VC_new;
        
      }
      
      epred_stan = dv_epred[i_obs];
      
      for(i in 1:n_obs){
        real epred_tmp = epred_stan[i];
        real sigma_tmp = sqrt(square(epred_tmp) * sigma_sq_p + sigma_sq_a + 
                              2*epred_tmp*sigma_p_a);
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

