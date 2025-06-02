// First Order Absorption (oral/subcutaneous)
// Two-compartment PK Model
// IIV on CL, VC, Q, and VP (full covariance matrix)
// proportional error - DV = IPRED*(1 + eps_p)
// Matrix-exponential solution
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   truncates the likelihood below at 0
// n_sim has to be extremely high (~ 5000-10000 for this model) for PSIS-LOO to
//   be useful
// Covariates: 
//   1) Body Weight on CL, VC, Q, and VP - (wt/70)^theta
//   2) Sex on VC (0 = male, 1 = female) - exp(theta*sex)
//   3) eGFR on CL (continuous) - (eGFR/90)^theta
//   4) Race on VC - exp(theta_vc_{race}*I(race == {race}))

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
  
  int<lower = 1> n_sim;
 
}
transformed data{ 
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 4; // Number of random effects
  int n_cmt = 2;    // Number of compartments - central, peripheral
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt);
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  real log_n_sim = log(n_sim); // This only needs calculated once for efficiency
                               // Probably won't make much difference
                               
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
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{
 
  vector[n_subjects] log_lik_subj;
 
  {
    
    vector[n_random] typical_values = to_vector({TVCL, TVVC, TVQ, TVVP});
    
    array[n_subjects] row_vector[n_sim] log_lik_subj_sim = 
      rep_array(rep_row_vector(0, n_sim), n_subjects);
    
    vector[n_races] theta_vc_race = [0, theta_vc_race2, theta_vc_race3, 
                                     theta_vc_race4]';
                                     
    vector[n_subjects] cl_p;
    vector[n_subjects] vc_p;
    vector[n_subjects] q_p;
    vector[n_subjects] vp_p;
    
    for(j in 1:n_subjects){
      
      real wt_adjustment_cl = wt_over_70[j]^theta_cl_wt;
      real wt_adjustment_vc = wt_over_70[j]^theta_vc_wt;
      real wt_adjustment_q = wt_over_70[j]^theta_q_wt;
      real wt_adjustment_vp = wt_over_70[j]^theta_vp_wt;
      real sex_adjustment_vc = exp(theta_vc_sex*sex[j]);
      real egfr_adjustment_cl = (egfr_over_90[j])^theta_cl_egfr;
      real race_adjustment_vc = exp(theta_vc_race[race[j]]);
      
      cl_p[j] = TVCL * wt_adjustment_cl * egfr_adjustment_cl;
      vc_p[j] = TVVC * wt_adjustment_vc * race_adjustment_vc * sex_adjustment_vc;
      q_p[j] = TVQ * wt_adjustment_q;
      vp_p[j] = TVVP * wt_adjustment_vp;
      
    }
    
    for(m in 1:n_sim){
      
      matrix[n_total, n_cmt] x_epred;
      vector[n_total] dv_epred;
      vector[n_obs] epred_stan;
      
      vector[n_random] eta_new = 
            multi_normal_cholesky_rng(rep_vector(0, n_random),
                                      diag_pre_multiply(omega, L)); 
      
      for(j in 1:n_subjects){
        
        real CL_new = cl_p[j]*exp(eta_new[1]);
        real VC_new = vc_p[j]*exp(eta_new[2]);
        real Q_new = q_p[j]*exp(eta_new[3]);
        real VP_new = vp_p[j]*exp(eta_new[4]);
        
        matrix[n_cmt, n_cmt] K_epred = rep_matrix(0, n_cmt, n_cmt);
        
        real ke_new = CL_new/VC_new;
        real k_cp_new = Q_new/VC_new;
        real k_pc_new = Q_new/VP_new;
    
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
        
                           
        dv_epred[subj_start[j]:subj_end[j]] = 
          x_epred[subj_start[j]:subj_end[j], 1] ./ VC_new;
        
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

