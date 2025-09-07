// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model with Indirect Response 3 (stimulation of production 
//   of response) PD Model
// IIV on CL, VC, and KA (full covariance matrix)
// IIV on KIN, KOUT, SC50, and SMAX (full covariance matrix) for PD. HILL is 
//   fixed to be 1
// proportional plus additive error for PK - DV = IPRED*(1 + eps_p_pk) + eps_a_pk
// proportional plus additive error for PD - DV = IPRED*(1 + eps_p_pd) + eps_a_pd
// General or coupled ODE solution using Torsten. The user can choose to get out
//    individual estimates of AUC, Cmax, and Tmax for PK and Tmax and Rmax for PD
// Predictions are generated from a normal that is truncated below at 0

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  vector depot_1cmt_ir3_ode(real t, vector y, array[] real params, 
                            array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real kin = params[4];
    real kout = params[5];
    real sc50 = params[6];
    real smax = params[7];   
    real hill = params[8];   // It's fixed to 1 in this particular model
    real r_0 = params[9];
    
    real ke = cl/vc;
    
    real conc = y[2]/vc;
    
    real stim = (smax*pow(conc, hill))/(pow(sc50, hill) + pow(conc, hill));
    real response = y[3] + r_0;
    
    vector[3] dydt;

    dydt[1] = -ka*y[1];
    dydt[2] = ka*y[1] - ke*y[2];
    dydt[3] = kin*(1 + stim) - kout*response;
    
    return dydt;
  }
  
  vector depot_1cmt_ir3_ode_coupled(real t, vector y, vector y_pk, 
                                    array[] real params, array[] real x_r, 
                                    array[] int x_i){
    
    real vc = params[2];

    real kin = params[4];
    real kout = params[5];
    real sc50 = params[6];
    real smax = params[7];   
    real hill = params[8];   // It's fixed to 1 in this particular model
    real r_0 = params[9];
    
    real conc = y_pk[2]/vc;
    
    real stim = (smax*pow(conc, hill))/(pow(sc50, hill) + pow(conc, hill));
    real response = y[1] + r_0;
    
    vector[1] dydt;

    dydt[1] = kin*(1 + stim) - kout*response;
    
    return dydt;
    
  }
  

  vector depot_1cmt_ir3_with_auc_cmax_ode(real t, vector y, array[] real params, 
                                          array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real kin = params[4];
    real kout = params[5];
    real sc50 = params[6];
    real smax = params[7];   // It's fixed to 1 in this particular model
    real hill = params[8];   // It's fixed to 1 in this particular model
    real r_0 = params[9];
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real ke = cl/vc;
    
    real conc = y[2]/vc;
    
    real stim = (smax*pow(conc, hill))/(pow(sc50, hill) + pow(conc, hill));
    real response = y[3] + r_0;
    
    real slope_pk = ka*y[1] - ke*y[2];
    real x_pk = slope_pk > 0 && conc > y[6] ? slope_pk/vc : 0;
    real z_pk = t <= t_1 || (slope_pk > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    real slope_pd = kin*(1 + stim) - kout*response;
    real x_pd = slope_pd > 0 && y[3] > y[8] ? slope_pd : 0;
    real z_pd = t <= t_1 || (slope_pd > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    vector[9] dydt;

    dydt[1] = -ka*y[1];
    dydt[2] = slope_pk;
    dydt[3] = slope_pd;
    dydt[4] = y[2];                                // AUC
    dydt[5] = t >= t_1 && t <= t_2 ? y[2] : 0;     // AUC_t_1-t_2
    dydt[6] = x_pk;                                // C_max 
    dydt[7] = z_pk;                                // t_max for PK
    dydt[8] = x_pd;                                // R_max
    dydt[9] = z_pd;                                // t_max for PD
    
    return dydt;
    
  }
  
}
data{
  
  int n_subjects;
  int n_subjects_new;
  int n_time_new;
  array[n_time_new] real time;
  array[n_time_new] real amt;
  array[n_time_new] int cmt;
  array[n_time_new] int evid;
  array[n_time_new] real rate;
  array[n_time_new] real ii;
  array[n_time_new] int addl;
  array[n_time_new] int ss;
  array[n_subjects_new] int subj_start;
  array[n_subjects_new] int subj_end;
  
  real<lower = 0> t_1;   // Time at which to start SS calculations (AUC_ss, C_max_ss, ...)
  real<lower = t_1> t_2; // Time at which to end SS calculations (AUC_ss, C_max_ss, ...)
  
  int<lower = 0, upper = 1> want_auc_cmax; // Want AUC, Cmax, Tmax, Rmin, Tmin?
  
}
transformed data{ 
  
  int n_random_pk = 3;
  int n_random_pd = 4;
  
  int n_cmt_pk = 2;    // number of ODEs in PK model (depot, central)
  int n_cmt_pd = 1;    // number of ODEs in PD system
  int n_cmt_extra = 6; // number of ODEs for AUC, Tmax, Rmin, ...
  int n_cmt = want_auc_cmax ? n_cmt_pk + n_cmt_pd + n_cmt_extra : n_cmt_pk + n_cmt_pd;
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);  // Hardcoding, but could be data or a parameter in another situation
  
  array[1, 2] real x_r = {{t_1, t_2}};
  
  real TVHILL = 1.0;                                                   // HILL is fixed to 1.0 in this model,
  vector<lower = 0>[n_subjects] HILL = rep_vector(TVHILL, n_subjects); // but it could be data or a parameter in another model.
                                                                       // Putting this here will require the least amount of 
                                                                       // changes in the code if that change is made

}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = TVCL/TVVC> TVKA;
  
  vector<lower = 0>[n_random_pk] omega_pk;
  cholesky_factor_corr[n_random_pk] L_pk;
  
  vector<lower = 0>[2] sigma_pk;
  cholesky_factor_corr[2] L_Sigma_pk;
  
  matrix[n_random_pk, n_subjects] Z_pk;
  
  real<lower = 0> TVKIN;       
  real<lower = 0> TVKOUT; 
  real<lower = 0> TVSC50;
  real<lower = 0> TVSMAX;
  
  vector<lower = 0>[n_random_pd] omega_pd;
  cholesky_factor_corr[n_random_pd] L_pd;
  
  vector<lower = 0>[2] sigma_pd;
  cholesky_factor_corr[2] L_Sigma_pd;
  
  matrix[n_random_pd, n_subjects] Z_pd;
  
}
generated quantities{

  vector[n_time_new] epred_stan; // f(TVs, x, eta = eta_new), eta_new ~ multi_normal(0, Omega) 
  vector[n_time_new] epred;      // epred_stan + error
  
  vector[want_auc_cmax ? n_time_new : 0] auc;          // AUC from 0 up until t
  vector[want_auc_cmax ? n_subjects_new : 0] auc_ss;   // AUC from t1 up to t2 (AUC_ss)
  vector[want_auc_cmax ? n_subjects_new : 0] c_max;    // Cmax between t1 and t2 (c_max_ss)
  vector[want_auc_cmax ? n_subjects_new : 0] t_max;    // Tmax between t1 and t2, then subtract off t1
  vector[want_auc_cmax ? n_subjects_new : 0] r_max;    // Maximum response
  vector[want_auc_cmax ? n_subjects_new : 0] t_max_pd; // Tmax for PD

  vector[n_subjects_new] CL;
  vector[n_subjects_new] VC;
  vector[n_subjects_new] KA;
  
  vector[n_subjects_new] KIN;
  vector[n_subjects_new] KOUT;
  vector[n_subjects_new] SC50;
  vector[n_subjects_new] SMAX;
  vector[n_subjects_new] R_0;
  
  vector[n_subjects_new] t_half;  

  {
    row_vector[n_random_pk] typical_values_pk = to_row_vector({TVCL, TVVC, TVKA});
    
    matrix[n_subjects_new, n_random_pk] eta_new_pk;
    matrix[n_subjects_new, n_random_pk] theta_new_pk;
    
    row_vector[n_random_pd] typical_values_pd = to_row_vector({TVKIN, TVKOUT, 
                                                               TVSC50, TVSMAX});
    
    matrix[n_subjects_new, n_random_pd] eta_new_pd;   
    matrix[n_subjects_new, n_random_pd] theta_new_pd; 
    
    matrix[2, 2] R_Sigma_pk = multiply_lower_tri_self_transpose(L_Sigma_pk);
    matrix[2, 2] Sigma_pk = quad_form_diag(R_Sigma_pk, sigma_pk);
    
    real sigma_sq_p_pk = Sigma_pk[1, 1];
    real sigma_sq_a_pk = Sigma_pk[2, 2];
    real sigma_p_a_pk = Sigma_pk[1, 2];
    
    matrix[2, 2] R_Sigma_pd = multiply_lower_tri_self_transpose(L_Sigma_pd);
    matrix[2, 2] Sigma_pd = quad_form_diag(R_Sigma_pd, sigma_pd);
    
    real sigma_sq_p_pd = Sigma_pd[1, 1];
    real sigma_sq_a_pd = Sigma_pd[2, 2];
    real sigma_p_a_pd = Sigma_pd[1, 2];
    
    matrix[n_time_new, n_cmt] x_epred;
    
    for(i in 1:n_subjects_new){
      eta_new_pk[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random_pk),
                                            diag_pre_multiply(omega_pk, L_pk))';
      eta_new_pd[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random_pd),
                                            diag_pre_multiply(omega_pd, L_pd))';
    }
    theta_new_pk = (rep_matrix(typical_values_pk, n_subjects_new) .* exp(eta_new_pk));
    theta_new_pd = (rep_matrix(typical_values_pd, n_subjects_new) .* exp(eta_new_pd));
    
    for(j in 1:n_subjects_new){
      
      row_vector[n_random_pk] theta_j_new_pk = theta_new_pk[j]; // access the PK parameters for subject j's epred
      row_vector[n_random_pd] theta_j_new_pd = theta_new_pd[j]; // access the PD parameters for subject j's epred
      
      CL[j] = theta_j_new_pk[1];
      VC[j] = theta_j_new_pk[2];
      KA[j] = theta_j_new_pk[3];
      t_half[j] = log(2)/(CL[j]/VC[j]);
      
      KIN[j] = theta_j_new_pd[1];
      KOUT[j] = theta_j_new_pd[2];
      SC50[j] = theta_j_new_pd[3];
      SMAX[j] = theta_j_new_pd[4];
      // HILL[j] = ...;
      R_0[j] = KIN[j]/KOUT[j];
      
      if(want_auc_cmax == 1){
        
        x_epred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(depot_1cmt_ir3_with_auc_cmax_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j], KA[j], 
                          KIN[j], KOUT[j], SC50[j], SMAX[j], HILL[j], R_0[j]},
                         bioav, tlag, x_r)';
                        
          
        auc[subj_start[j]:subj_end[j]] = 
                                x_epred[subj_start[j]:subj_end[j], 4] ./ VC[j];
          
        auc_ss[j] = max(x_epred[subj_start[j]:subj_end[j], 5]) / VC[j];
        c_max[j] = max(x_epred[subj_start[j]:subj_end[j], 6]);
        t_max[j] = max(x_epred[subj_start[j]:subj_end[j], 7]) - t_1;
        
        r_max[j] = max(x_epred[subj_start[j]:subj_end[j], 8]) + R_0[j];
        t_max_pd[j] = max(x_epred[subj_start[j]:subj_end[j], 9]) - t_1;
        
      }else{
        
        x_epred[subj_start[j]:subj_end[j],] =
          pmx_solve_onecpt_rk45(depot_1cmt_ir3_ode_coupled,
                                n_cmt_pd,
                                time[subj_start[j]:subj_end[j]],
                                amt[subj_start[j]:subj_end[j]],
                                rate[subj_start[j]:subj_end[j]],
                                ii[subj_start[j]:subj_end[j]],
                                evid[subj_start[j]:subj_end[j]],
                                cmt[subj_start[j]:subj_end[j]],
                                addl[subj_start[j]:subj_end[j]],
                                ss[subj_start[j]:subj_end[j]],
                                {CL[j], VC[j], KA[j],
                                 KIN[j], KOUT[j], SC50[j],
                                 SMAX[j], HILL[j], R_0[j]},
                                bioav, tlag)';
        
      }
      
      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){
          epred_stan[k] = x_epred[k, 2] / VC[j];
        }else if(cmt[k] == 3){
          epred_stan[k] = x_epred[k, 3] + R_0[j];
        }
      }
    
    }
    

    for(i in 1:n_time_new){
      if(epred_stan[i] == 0){
        epred[i] = 0;
      }else{
        if(cmt[i] == 2 || cmt[i] == 3){
          real epred_tmp = epred_stan[i];
          real sigma_tmp_e = cmt[i] == 2 ? 
            sqrt(square(epred_tmp) * sigma_sq_p_pk + sigma_sq_a_pk + 
                              2*epred_tmp*sigma_p_a_pk) : 
            sqrt(square(epred_tmp) * sigma_sq_p_pd + sigma_sq_a_pd + 
                              2*epred_tmp*sigma_p_a_pd);
          epred[i] = normal_lb_rng(epred_tmp, sigma_tmp_e, 0.0);
        }
      }
    }
  
  }

}


