// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model, Friberg-Karlsson neutropenia PD model
// IIV on CL, VC, and Ka (full covariance matrix)
// IIV on MTT, CIRC0, GAMMA, ALPHA (full covariance matrix) for PD
// exponential error on PK - DV = IPRED*exp(eps)
// exponential error on PD - DV = IPRED*exp(eps_pd)
// General ODE solution using Torsten to get out individual estimates of AUC, 
//   Cmax, Tmax, ...
// Predictions are generated from a normal that is truncated below at 0

functions{
  
  vector depot_1cmt_friberg_ode(real t, vector y, array[] real params, 
                                array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real mtt = params[4];
    real circ_0 = params[5];
    real gamma = params[6];
    real alpha = params[7];
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real ke = cl/vc;
    real k_tr = 4/mtt; // k_tr = (n_tr + 1)/mtt    
    
    real conc = y[2]/vc;
    
    real e_drug = fmin(1.0, alpha*conc); // Maybe reparameterize this so no more fmin?
    real prol = y[3] + circ_0;
    real transit_1 = y[4] + circ_0; 
    real transit_2 = y[5] + circ_0;
    real transit_3 = y[6] + circ_0;
    // real circ = y[7] + circ_0; // fmax(machine_precision(), y[7] + circ_0)
    real circ = fmax(machine_precision(), y[7] + circ_0);
    
    real slope_pk = ka*y[1] - ke*y[2];
    real x_pk = slope_pk > 0 && conc > y[10] ? slope_pk/vc : 0;
    real z_pk = t <= t_1 || (slope_pk > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    real slope_pd = k_tr*(transit_3 - circ);
    real x_pd = slope_pd < 0 && y[7] < y[12] ? slope_pd : 0;
    
    vector[12] dydt;
    
    dydt[1] = -ka*y[1];                               // depot
    dydt[2] = slope_pk;                               // central
    dydt[3] = k_tr*prol*((1 - e_drug)*(circ_0/circ)^gamma - 1);  // proliferative cells
    dydt[4] = k_tr*(prol - transit_1);                // transit 1
    dydt[5] = k_tr*(transit_1 - transit_2);           // transit 2
    dydt[6] = k_tr*(transit_2 - transit_3);           // transit 3
    dydt[7] = slope_pd;                               // circulating blood cells
    dydt[8] = y[2];                                   // AUC
    dydt[9] = t >= t_1 && t <= t_2 ? y[2] : 0;        // AUC_t_1-t_2
    dydt[10] = x_pk;                                  // C_max 
    dydt[11] = z_pk;                                  // t_max for PK
    dydt[12] = x_pd;                                  // R_min
    
    return dydt;
  }
  
}
data{
  
  int n_subjects;
  int n_time_new;
  array[n_time_new] real time;
  array[n_time_new] real amt;
  array[n_time_new] int cmt;
  array[n_time_new] int evid;
  array[n_time_new] real rate;
  array[n_time_new] real ii;
  array[n_time_new] int addl;
  array[n_time_new] int ss;
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;
  
  real<lower = 0> t_1;   // Time at which to start SS calculations (AUC_ss, C_max_ss, ...)
  real<lower = t_1> t_2; // Time at which to end SS calculations (AUC_ss, C_max_ss, ...)
 
}
transformed data{ 
  
  int n_random = 3;
  int n_random_pd = 4;
  
  int n_cmt = 2;       // number of ODEs in PK model (depot, central)
  int n_cmt_pd = 5;    // number of ODEs in PD system
  int n_cmt_extra = 5; // number of ODEs for AUC, Tmax, Rmin, ...
  
  array[n_cmt + n_cmt_pd + n_cmt_extra] real bioav = 
                                 rep_array(1.0, n_cmt + n_cmt_pd + n_cmt_extra); 
  array[n_cmt + n_cmt_pd + n_cmt_extra] real tlag = 
                                 rep_array(0.0, n_cmt + n_cmt_pd + n_cmt_extra);
  
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
  
  real<lower = 0> TVMTT;       
  real<lower = 0> TVCIRC0; 
  real<lower = 0> TVGAMMA;
  real<lower = 0> TVALPHA;
  
  vector<lower = 0>[n_random_pd] omega_pd;
  cholesky_factor_corr[n_random_pd] L_pd;
  
  real<lower = 0> sigma_pd;
  
  matrix[n_random_pd, n_subjects] Z_pd;
  
}
generated quantities{
  
  vector[n_time_new] ipred;   // ipred for the observed individuals at the new timepoints
  vector[n_time_new] pred;    // pred for the observed individuals at the new timepoints
  vector[n_time_new] dv;      // dv for the observed individuals at the new timepoints
  vector[n_time_new] auc;     // auc for the observed individuals from time 0 to the new timepoint
  vector[n_subjects] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[n_subjects] c_max;   // Cmax between t1 and t2 (c_max_ss)
  vector[n_subjects] t_max;   // Tmax between t1 and t2, then subtract off t1
  vector[n_subjects] t_half;  // half-life
  vector[n_subjects] r_min;   // Minimum response
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  
  vector[n_subjects] MTT;
  vector[n_subjects] CIRC0;
  vector[n_subjects] GAMMA;
  vector[n_subjects] ALPHA;
 
  {
    
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA});

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    row_vector[n_random_pd] typical_values_pd = to_row_vector({TVMTT, TVCIRC0, 
                                                               TVGAMMA, TVALPHA});
    
    matrix[n_subjects, n_random_pd] eta_pd = 
                                      diag_pre_multiply(omega_pd, L_pd * Z_pd)';  
    matrix[n_subjects, n_random_pd] theta_pd =
                     (rep_matrix(typical_values_pd, n_subjects) .* exp(eta_pd));
    
    matrix[n_time_new, n_cmt + n_cmt_pd + n_cmt_extra] x_pred;
    matrix[n_time_new, n_cmt + n_cmt_pd + n_cmt_extra] x_ipred;
    
    CL = col(theta, 1);
    VC = col(theta, 2);
    KA = col(theta, 3);
    
    MTT = col(theta_pd, 1);
    CIRC0 = col(theta_pd, 2);
    GAMMA = col(theta_pd, 3);
    ALPHA = col(theta_pd, 4);
    
    for(j in 1:n_subjects){
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_1cmt_friberg_ode,
                       n_cmt + n_cmt_pd + n_cmt_extra,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {CL[j], VC[j], KA[j], 
                          MTT[j], CIRC0[j], GAMMA[j], ALPHA[j]},
                       bioav, tlag, x_r)';

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_1cmt_friberg_ode,
                       n_cmt + n_cmt_pd + n_cmt_extra,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {TVCL, TVVC, TVKA, 
                        TVMTT, TVCIRC0, TVGAMMA, TVALPHA},
                       bioav, tlag, x_r)';

      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){
          // ipred[k] = x_ipred[k, 2] / VC[j];
          // pred[k] = x_pred[k, 2] / TVVC;
          ipred[k] = fmax(machine_precision(), x_ipred[k, 2] / VC[j]);
          pred[k] = fmax(machine_precision(), x_pred[k, 2] / TVVC);
        }else if(cmt[k] == 4){
          // ipred[k] = x_ipred[k, 7] + CIRC0[j];
          // pred[k] = x_pred[k, 7] + TVCIRC0;
          ipred[k] = fmax(machine_precision(), x_ipred[k, 7] + CIRC0[j]);
          pred[k] = fmax(machine_precision(), x_pred[k, 7] + TVCIRC0);
        }
      }
      
      auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 8] ./ VC[j];  
      auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], 9]) / VC[j];
      c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 10]);
      t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 11]) - t_1;
      t_half[j] = log(2)/(CL[j]/VC[j]);
      
      r_min[j] = min(x_ipred[subj_start[j]:subj_end[j], 12]) + CIRC0[j];
    
    }
    
    for(i in 1:n_time_new){
      if(ipred[i] == 0){
        dv[i] = 0;
      }else if(cmt[i] == 2 || cmt[i] == 4){
        real log_ipred_tmp = log(ipred[i]);
        real sigma_tmp = cmt[i] == 2 ? sigma : sigma_pd;
    
        dv[i] = lognormal_rng(log_ipred_tmp, sigma_tmp);
      }
    }
  }
}
