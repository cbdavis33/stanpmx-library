// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model with linear and Michaelis-Menten elimination
// IIV on CL, VC, VMAX, KM, KA (full covariance matrix)
// exponential error - DV = IPRED*exp(eps)
// General ODE solution using Torsten. The user can choose to get out individual
//    estimates of AUC, Cmax, and Tmax

functions{

  vector depot_1cmt_linear_plus_mm_ode(real t, vector y, array[] real params, 
                                       array[] real x_r, array[] int x_i){
  
    real cl = params[1];
    real vc = params[2];
    real vmax = params[3];
    real km = params[4];
    real ka = params[5];
    real ke = cl/vc;
    
    real conc = y[2]/vc;
    
    vector[2] dydt;

    dydt[1] = -ka*y[1];
    dydt[2] = ka*y[1] - ke*y[2] - vmax*conc/(km + conc);
    
    return dydt;
  }
  
  vector depot_1cmt_linear_plus_mm_with_auc_ode(real t, vector y, array[] real params, 
                                                array[] real x_r, array[] int x_i){
  
    real cl = params[1];
    real vc = params[2];
    real vmax = params[3];
    real km = params[4];
    real ka = params[5];
    real ke = cl/vc;
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real conc = y[2]/vc;
    
    real slope = ka*y[1] - ke*y[2] - vmax*conc/(km + conc);
    real x = slope > 0 && conc > y[5] ? slope/vc : 0;
    real z = t <= t_1 || (slope > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    vector[6] dydt;
    
    dydt[1] = -ka*y[1];
    dydt[2] = slope;
    dydt[3] = y[2];                                // AUC
    dydt[4] = t >= t_1 && t <= t_2 ? y[2] : 0;     // AUC_t_1-t_2
    dydt[5] = x;                                   // C_max
    dydt[6] = z;                                   // t_max
    
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
  
  int<lower = 0, upper = 1> want_auc_cmax; // Want AUC, Cmax, and Tmax? 
 
}
transformed data{ 
  
  int n_random = 5;
  int n_cmt = want_auc_cmax ? 6 : 2; // Number of compartments - depot, central (AUC, AUC_t1_t2, Cmax_ss, Tmax_ss))
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt);
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  array[1, 2] real x_r = {{t_1, t_2}};

}
parameters{ 
  
  real<lower = 0> TVCL;
  real<lower = 0> TVVC;
  real<lower = 0> TVVMAX;       
  real<lower = 0> TVKM;
  real<lower = 0> TVKA;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{
  
  vector[n_time_new] pred;       // f(TVs, x, eta = 0) 
  vector[n_time_new] epred_stan; // f(TVs, x, eta = eta_new), eta_new ~ multi_normal(0, Omega) 
  vector[n_time_new] epred;      // epred_stan + error
  vector[n_time_new] ipred;      // f(TVs, x, eta = eta_i), eta_i are etas for observed subjects
  vector[n_time_new] dv;         // ipred + error
  
  vector[want_auc_cmax ? n_time_new : 0] auc;     // AUC from 0 up until t
  vector[want_auc_cmax ? n_subjects : 0] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[want_auc_cmax ? n_subjects : 0] c_max;   // Cmax between t1 and t2 (c_max_ss)
  vector[want_auc_cmax ? n_subjects : 0] t_max;   // Tmax between t1 and t2, then subtract off t1
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] VMAX;
  vector[n_subjects] KM;
  vector[n_subjects] KA;
 
  {
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVVMAX, TVKM, TVKA});

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    matrix[n_subjects, n_random] eta_new;
    matrix[n_subjects, n_random] theta_new;
    
    for(i in 1:n_subjects){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects) .* exp(eta_new));
    
    vector[n_subjects] CL_new;
    vector[n_subjects] VC_new;
    vector[n_subjects] VMAX_new;
    vector[n_subjects] KM_new;
    vector[n_subjects] KA_new;

    matrix[n_time_new, 2] x_pred;
    matrix[n_time_new, 2] x_epred;
    matrix[n_time_new, n_cmt] x_ipred;
    
    for(j in 1:n_subjects){
      
      row_vector[n_random] theta_j = theta[j]; // access the parameters for subject j
      row_vector[n_random] theta_j_new = theta_new[j]; // access the parameters for subject j's epred
      
      real cl_p = TVCL;
      real vc_p = TVVC;
      real vmax_p = TVVMAX;
      real km_p = TVKM;
      real ka_p = TVKA;
      
      CL[j] = theta_j[1];
      VC[j] = theta_j[2];
      VMAX[j] = theta_j[3];
      KM[j] = theta_j[4];
      KA[j] = theta_j[5];
      
      CL_new[j] = theta_j_new[1];
      VC_new[j] = theta_j_new[2];
      VMAX_new[j] = theta_j_new[3];
      KM_new[j] = theta_j_new[4];
      KA_new[j] = theta_j_new[5];
      
      x_pred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(depot_1cmt_linear_plus_mm_ode,
                         2,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {cl_p, vc_p, vmax_p, km_p, ka_p}, 
                         bioav[1:2], tlag[1:2], x_r)';

      pred[subj_start[j]:subj_end[j]] = 
        x_pred[subj_start[j]:subj_end[j], 2] ./ vc_p;
      
      x_epred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(depot_1cmt_linear_plus_mm_ode,
                         2,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL_new[j], VC_new[j], VMAX_new[j], KM_new[j], KA_new[j]}, 
                         bioav[1:2], tlag[1:2], x_r)';
      
      epred_stan[subj_start[j]:subj_end[j]] = 
        x_epred[subj_start[j]:subj_end[j], 2] ./ VC_new[j];
      
      
      if(want_auc_cmax == 1){
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(depot_1cmt_linear_plus_mm_with_auc_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j], VMAX[j], KM[j], KA[j]}, 
                         bioav, tlag, x_r)';
                        
        ipred[subj_start[j]:subj_end[j]] = 
          x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
          
        auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 3] ./ VC[j];
          
        auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], 4]) / VC[j];
        c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]);
        t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 6]) - t_1;
          
      }else{
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(depot_1cmt_linear_plus_mm_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j], VMAX[j], KM[j], KA[j]}, 
                         bioav, tlag, x_r)';

        ipred[subj_start[j]:subj_end[j]] = 
          x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
        
      }

    }
  
    for(i in 1:n_time_new){
      
      if(epred_stan[i] == 0){
        epred[i] = 0;
      }else{
        epred[i] = lognormal_rng(log(epred_stan[i]), sigma);
      }
      
      if(ipred[i] == 0){
        dv[i] = 0;
      }else{
        dv[i] = lognormal_rng(log(ipred[i]), sigma);
      } 
      
    }
  }
}
