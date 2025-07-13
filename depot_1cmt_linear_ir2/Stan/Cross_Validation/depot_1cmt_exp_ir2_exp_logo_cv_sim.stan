// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model with Indirect Response 2 (inhibition of dissipation 
//   of response) PD Model
// IIV on CL, VC, and KA (full covariance matrix)
// IIV on KIN, KOUT, IC50 (full covariance matrix) for PD. IMAX and HILL are 
//   fixed to be 1
// exponential error on PK - DV = IPRED*exp(eps_pk)
// exponential error on PD - DV = IPRED*exp(eps_pd)
// n_sim has to be extremely high (~ 5000-10000 for this model) for PSIS-LOO to
//   be useful

functions{
  
  vector depot_1cmt_ir2_ode(real t, vector y, array[] real params, 
                            array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real kin = params[4];
    real kout = params[5];
    real ic50 = params[6];
    real imax = params[7];   // It's fixed to 1 in this particular model
    real hill = params[8];   // It's fixed to 1 in this particular model
    real r_0 = params[9];
    
    real ke = cl/vc;
    
    real conc = y[2]/vc;
    
    real inh = (imax*pow(conc, hill))/(pow(ic50, hill) + pow(conc, hill));
    real response = y[3] + r_0;
    
    vector[3] dydt;

    dydt[1] = -ka*y[1];
    dydt[2] = ka*y[1] - ke*y[2];
    dydt[3] = kin - kout*(1 - inh)*response;
    
    return dydt;
    
  }
  
  vector depot_1cmt_ir2_ode_coupled(real t, vector y, vector y_pk, 
                                    array[] real params, array[] real x_r, 
                                    array[] int x_i){
    
    real vc = params[2];

    real kin = params[4];
    real kout = params[5];
    real ic50 = params[6];
    real imax = params[7];   // It's fixed to 1 in this particular model
    real hill = params[8];   // It's fixed to 1 in this particular model
    real r_0 = params[9];
    
    real conc = y_pk[2]/vc;
    
    real inh = (imax*pow(conc, hill))/(pow(ic50, hill) + pow(conc, hill));
    real response = y[1] + r_0;
    
    vector[1] dydt;

    dydt[1] = kin - kout*(1 - inh)*response;
    
    return dydt;
    
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
  
  int<lower = 1> n_sim;
 
}
transformed data{ 
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random_pk = 3;
  int n_random_pd = 3;
  
  int n_cmt_pk = 2;  // number of ODEs in PK model (depot, central)
  int n_cmt_pd = 1;  // number of ODEs in PD system
  int n_cmt = n_cmt_pk + n_cmt_pd;
  
  array[n_cmt_pk + n_cmt_pd] real bioav = rep_array(1.0, n_cmt_pk + n_cmt_pd); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt_pk + n_cmt_pd] real tlag = rep_array(0.0, n_cmt_pk + n_cmt_pd);
  
  real log_n_sim = log(n_sim); // This only needs calculated once for efficiency
                               // Probably won't make much difference
                               
  real TVIMAX = 1.0;
  real TVHILL = 1.0;
  vector<lower = 0>[n_subjects] IMAX = rep_vector(TVIMAX, n_subjects); // IMAX and HILL are both fixed to 1.0 in this model,
  vector<lower = 0>[n_subjects] HILL = rep_vector(TVHILL, n_subjects); // but thry could be data or a parameter in another model.
                                                                       // Putting this here will require the least amount of 
                                                                       // changes in the code if that change is made

}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = TVCL/TVVC> TVKA;
  
  vector<lower = 0>[n_random_pk] omega_pk;
  cholesky_factor_corr[n_random_pk] L_pk;
  
  real<lower = 0> sigma_pk;
  
  matrix[n_random_pk, n_subjects] Z_pk;
  
  real<lower = 0> TVKIN;       
  real<lower = 0> TVKOUT; 
  real<lower = 0> TVIC50;
  
  vector<lower = 0>[n_random_pd] omega_pd;
  cholesky_factor_corr[n_random_pd] L_pd;
  
  real<lower = 0> sigma_pd;
  
  matrix[n_random_pd, n_subjects] Z_pd;
  
}
generated quantities{
 
  vector[n_subjects] log_lik_subj;
 
  {
    
    vector[n_random_pk] typical_values_pk = to_vector({TVCL, TVVC, TVKA});
    vector[n_random_pd] typical_values_pd = to_vector({TVKIN, TVKOUT, TVIC50});
    array[n_subjects] row_vector[n_sim] log_lik_subj_sim = 
      rep_array(rep_row_vector(0, n_sim), n_subjects);
    
    for(m in 1:n_sim){
      
      matrix[n_total, n_cmt] x_epred;
      vector[n_total] dv_epred;
      vector[n_obs] epred_stan;
      
      vector[n_random_pk] eta_new_pk = 
            multi_normal_cholesky_rng(rep_vector(0, n_random_pk),
                                      diag_pre_multiply(omega_pk, L_pk)); 
      vector[n_random_pk] theta_new_pk = typical_values_pk .* exp(eta_new_pk);
      
      vector[n_random_pd] eta_new_pd = 
            multi_normal_cholesky_rng(rep_vector(0, n_random_pd),
                                      diag_pre_multiply(omega_pd, L_pd)); 
      vector[n_random_pd] theta_new_pd = typical_values_pd .* exp(eta_new_pd);
      
      for(j in 1:n_subjects){
        
        // // placeholder for making adjustments for covariates if necessary
        // real wt_over_70 = wt[j]/70;
        // real wt_adjustment_cl = wt_over_70^theta_cl_wt;
        // real wt_adjustment_vc = wt_over_70^theta_vc_wt;
        
        real CL_new = theta_new_pk[1]; // * weight_adjustment_cl;
        real VC_new = theta_new_pk[2]; // * weight_adjustment_vc;
        real KA_new = theta_new_pk[3];
        
        real KIN_new = theta_new_pd[1];  
        real KOUT_new = theta_new_pd[2]; 
        real IC50_new = theta_new_pd[3];
        real IMAX_new = IMAX[j];
        real HILL_new = HILL[j];
        real R_0_new = KIN_new/KOUT_new;
        
        x_epred[subj_start[j]:subj_end[j],] =
          pmx_solve_onecpt_rk45(depot_1cmt_ir2_ode_coupled,
                                n_cmt_pd,
                                time[subj_start[j]:subj_end[j]],
                                amt[subj_start[j]:subj_end[j]],
                                rate[subj_start[j]:subj_end[j]],
                                ii[subj_start[j]:subj_end[j]],
                                evid[subj_start[j]:subj_end[j]],
                                cmt[subj_start[j]:subj_end[j]],
                                addl[subj_start[j]:subj_end[j]],
                                ss[subj_start[j]:subj_end[j]],
                                {CL_new, VC_new, KA_new,
                                 KIN_new, KOUT_new, IC50_new,
                                 IMAX_new, HILL_new, R_0_new},
                                bioav, tlag)';
        
        // x_epred[subj_start[j]:subj_end[j],] =
        //   pmx_solve_rk45(depot_1cmt_ir2_ode,
        //                  n_cmt_pk + n_cmt_pd,
        //                  time[subj_start[j]:subj_end[j]],
        //                  amt[subj_start[j]:subj_end[j]],
        //                  rate[subj_start[j]:subj_end[j]],
        //                  ii[subj_start[j]:subj_end[j]],
        //                  evid[subj_start[j]:subj_end[j]],
        //                  cmt[subj_start[j]:subj_end[j]],
        //                  addl[subj_start[j]:subj_end[j]],
        //                  ss[subj_start[j]:subj_end[j]],
        //                  {CL_new, VC_new, KA_new, 
        //                   KIN_new, KOUT_new, IC50_new, 
        //                   IMAX_new, HILL_new, R_0_new},
        //                  bioav, tlag)';
                           
        for(k in subj_start[j]:subj_end[j]){
          if(cmt[k] == 2){
            dv_epred[k] = x_epred[k, 2] / VC_new;
          }else if(cmt[k] == 3){
            dv_epred[k] = x_epred[k, 3] + R_0_new;
          }
        }
        
      }
      
      epred_stan = dv_epred[i_obs];
      
      for(i in 1:n_obs){
        
        if(cmt[i_obs[i]] == 2 || cmt[i_obs[i]] == 3){
          real sigma_tmp = cmt[i_obs[i]] == 2 ? sigma_pk : sigma_pd;
          if(bloq_obs[i] == 1){
            log_lik_subj_sim[dv_obs_id[i], m] += 
              lognormal_lcdf(lloq_obs[i] | log(epred_stan[i]), sigma_tmp);
          }else{
            log_lik_subj_sim[dv_obs_id[i], m] += 
              lognormal_lpdf(dv_obs[i] | log(epred_stan[i]), sigma_tmp);
          }

        }
     
      }
    }
    for(j in 1:n_subjects){
      log_lik_subj[j] = log_sum_exp(log_lik_subj_sim[j]) - log_n_sim;
    }
    
  }
}

