// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model with Indirect Response 2 PD Model
// Inhibitory Effect Compartment PD Model
// IIV on CL, VC, and KA (full covariance matrix)
// IIV on E0, KE0, and EC50 (full covariance matrix) for PD. EMAX and HILL are 
//   fixed to be 1
// For PD - IPRED = E0*(1 - emax*conc_effcmt^hill/(ec50^hill + conc_effcmt^hill))
// exponential error on PK - DV = IPRED*exp(eps)
// exponential error on PD - DV = IPRED*exp(eps_pd)
// General ODE solution using Torsten to get out individual estimates of AUC, 
//   Cmax, Tmax, ...
// Predictions are generated from a normal that is truncated below at 0

functions{
  
  vector depot_1cmt_effcmt_ode(real t, vector y, array[] real params, 
                               array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real ke0 = params[4];
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real ke = cl/vc;
    
    real conc = y[2]/vc;
    
    real slope_pk = ka*y[1] - ke*y[2];
    real x_pk = slope_pk > 0 && conc > y[6] ? slope_pk/vc : 0;
    real z_pk = t <= t_1 || (slope_pk > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    vector[7] dydt;

    dydt[1] = -ka*y[1];
    dydt[2] = slope_pk;
    dydt[3] = ke0*(y[2] - y[3]);
    dydt[4] = y[2];                                // AUC
    dydt[5] = t >= t_1 && t <= t_2 ? y[2] : 0;     // AUC_t_1-t_2
    dydt[6] = x_pk;                                // C_max 
    dydt[7] = z_pk;                                // t_max for PK
    
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
  
}
transformed data{ 
  
  int n_random = 3;
  int n_random_pd = 3;
  
  int n_cmt = 2;       // number of ODEs in PK model (depot, central)
  int n_cmt_pd = 1;    // number of ODEs in PD system
  int n_cmt_extra = 4; // number of ODEs for AUC, Tmax, Cmax, ...
  
  array[n_cmt + n_cmt_pd + n_cmt_extra] real bioav = 
                                 rep_array(1.0, n_cmt + n_cmt_pd + n_cmt_extra); 
  array[n_cmt + n_cmt_pd + n_cmt_extra] real tlag = 
                                 rep_array(0.0, n_cmt + n_cmt_pd + n_cmt_extra);
  
  real TVEMAX = 1.0;
  real TVHILL = 1.0;
  vector<lower = 0>[n_subjects] EMAX = rep_vector(TVEMAX, n_subjects); // EMAX and HILL are both fixed to 1.0 in this model,
  vector<lower = 0>[n_subjects] HILL = rep_vector(TVHILL, n_subjects); // but it could be data or a parameter in another model.
             
  array[1, 2] real x_r = {{t_1, t_2}};

}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = TVCL/TVVC> TVKA;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma;
  
  matrix[n_random, n_subjects] Z;
  
  real<lower = 0> TVE0;       
  real<lower = 0> TVKE0; 
  real<lower = 0> TVEC50;
  
  vector<lower = 0>[n_random_pd] omega_pd;
  cholesky_factor_corr[n_random_pd] L_pd;
  
  real<lower = 0> sigma_pd;
  
  matrix[n_random_pd, n_subjects] Z_pd;
  
}
generated quantities{

  vector[n_time_new] ipred;       // ipred for the observed individuals at the new timepoints
  vector[n_time_new] pred;        // pred for the observed individuals at the new timepoints
  vector[n_time_new] dv;          // dv for the observed individuals at the new timepoints
  vector[n_time_new] auc;         // auc for the observed individuals from time 0 to the new timepoint
  vector[n_subjects_new] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[n_subjects_new] c_max;   // Cmax between t1 and t2 (c_max_ss)
  vector[n_subjects_new] t_max;   // Tmax between t1 and t2, then subtract off t1
  vector[n_subjects_new] t_half;  // half-life
  
  vector[n_subjects_new] CL;
  vector[n_subjects_new] VC;
  vector[n_subjects_new] KA;
  
  vector[n_subjects_new] E0;
  vector[n_subjects_new] KE0;
  vector[n_subjects_new] EC50;

  {
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA});
    
    matrix[n_subjects_new, n_random] eta_new;
    matrix[n_subjects_new, n_random] theta_new;
    
    row_vector[n_random_pd] typical_values_pd = to_row_vector({TVE0, TVKE0, 
                                                               TVEC50});
    
    matrix[n_subjects_new, n_random_pd] eta_new_pd;   
    matrix[n_subjects_new, n_random_pd] theta_new_pd; 
    
    matrix[n_time_new, n_cmt + n_cmt_pd + n_cmt_extra] x_ipred;
    matrix[n_time_new, n_cmt + n_cmt_pd + n_cmt_extra] x_pred;

    for(i in 1:n_subjects_new){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
                                               
      eta_new_pd[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random_pd),
                                            diag_pre_multiply(omega_pd, L_pd))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects_new) .* exp(eta_new));
    theta_new_pd = 
      (rep_matrix(typical_values_pd, n_subjects_new) .* exp(eta_new_pd));
    
    CL = col(theta_new, 1);
    VC = col(theta_new, 2);
    KA = col(theta_new, 3);
    
    E0 = col(theta_new_pd, 1);
    KE0 = col(theta_new_pd, 2);
    EC50 = col(theta_new_pd, 3);
    
    for(j in 1:n_subjects_new){
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_1cmt_effcmt_ode,
                       n_cmt + n_cmt_pd + n_cmt_extra,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {CL[j], VC[j], KA[j], KE0[j]},
                       bioav, tlag, x_r)';

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_1cmt_effcmt_ode,
                       n_cmt + n_cmt_pd + n_cmt_extra,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {TVCL, TVVC, TVKA, TVKE0},
                       bioav, tlag, x_r)';

      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){
          ipred[k] = x_ipred[k, 2] / VC[j];
          pred[k] = x_pred[k, 2] / TVVC;
        }else if(cmt[k] == 3){
          real conc_effcmt = x_ipred[k, 3]/VC[j];
          real inh = (EMAX[j]*pow(conc_effcmt, HILL[j]))/(pow(EC50[j], HILL[j]) + 
                                                    pow(conc_effcmt, HILL[j]));
          ipred[k] = E0[j]*(1 - inh);
          
          real tv_conc_effcmt = x_pred[k, 3]/TVVC;
          real tvinh = (TVEMAX*pow(tv_conc_effcmt, TVHILL))/(pow(TVEC50, TVHILL) + 
                                                    pow(tv_conc_effcmt, TVHILL));
          pred[k] = TVE0*(1 - tvinh);
        }
      }
      
      auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 4] ./ VC[j];  
      auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]) / VC[j];
      c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 6]);
      t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 7]) - t_1;
      t_half[j] = log(2)/(CL[j]/VC[j]);
    
    }

    for(i in 1:n_time_new){
      if(ipred[i] == 0){
        dv[i] = 0;
      }else if(cmt[i] == 2 || cmt[i] == 3){
        real log_ipred_tmp = log(ipred[i]);
        real sigma_tmp = cmt[i] == 2 ? sigma : sigma_pd;
    
        dv[i] = lognormal_rng(log_ipred_tmp, sigma_tmp);
      }
    }
  
  }

}


