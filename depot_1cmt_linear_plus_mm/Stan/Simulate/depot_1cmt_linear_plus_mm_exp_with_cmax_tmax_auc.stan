// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model with Michaelis-Menten elimination
// IIV on VC, VMAX, KM, KA
// exponential error - DV = IPRED*exp(eps)
// ODE solution using Torsten

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
  int n_total;                  
  
  array[n_total] real amt;
  array[n_total] int cmt;
  array[n_total] int evid;
  array[n_total] real rate;
  array[n_total] real ii;
  array[n_total] int addl;
  array[n_total] int ss;
  array[n_total] real time;
  
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;
  
  real<lower = 0> TVCL;
  real<lower = 0> TVVC;
  real<lower = 0> TVVMAX;
  real<lower = 0> TVKM;
  real<lower = 0> TVKA;
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_vmax;
  real<lower = 0> omega_km;
  real<lower = 0> omega_ka;
  
  corr_matrix[5] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_vc_vmax, ... and then construct the 
                     // correlation matrix in transformed data, but it's easy
                     // enough to do in R
  
  real<lower = 0> sigma;
  
  int<lower = 1, upper = 2> solver; // 1 = rk45, 2 = bdf
  
  real<lower = 0> t_1;
  real<lower = t_1> t_2;
  
  int<lower = 0, upper = 1> want_auc_cmax; // Want AUC Cmax, and Tmax? 
                                           
}
transformed data{
  
  int n_random = 5;
  int n_cmt = want_auc_cmax ? 6 : 2; // Number of compartments - depot, central (AUC, AUC_t1_t2, Cmax_ss, Tmax_ss))

  vector[n_random] omega = [omega_cl, omega_vc, omega_vmax, omega_km, omega_ka]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);  // Hardcoding, but could be data or a parameter in another situation
  
  array[1, 2] real x_r = {{t_1, t_2}};
  
}
model{
  
}
generated quantities{
  
  vector[n_total] ipred; // concentration with no residual error
  vector[n_total] dv;    // concentration with residual error

  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] VMAX;
  vector[n_subjects] KM;
  vector[n_subjects] KA;
  
  vector[n_total] auc;            // AUC 
  vector[n_subjects] auc_t1_t2;   // AUC from t1 up to t2
  vector[n_subjects] c_max;       // Cmax
  vector[n_subjects] t_max;       // Tmax
  
  {
  
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVVMAX, TVKM, TVKA});
    
    matrix[n_subjects, n_random] eta;   
    matrix[n_subjects, n_random] theta; 
    
    matrix[n_total, n_cmt] x_ipred;
    
    for(i in 1:n_subjects){
      eta[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                           diag_pre_multiply(omega, L))';
    }
    theta = (rep_matrix(typical_values, n_subjects) .* exp(eta));

    CL = col(theta, 1);
    VC = col(theta, 2);
    VMAX = col(theta, 3);
    KM = col(theta, 4);
    KA = col(theta, 5);
    
    for(j in 1:n_subjects){
      
      if(want_auc_cmax){
        if(solver == 1){
          
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
                           
        }else{
          
          x_ipred[subj_start[j]:subj_end[j],] =
            pmx_solve_bdf(depot_1cmt_linear_plus_mm_with_auc_ode,
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
                         
        }
        
        auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 3] ./ VC[j];
      
        auc_t1_t2[j] = max(x_ipred[subj_start[j]:subj_end[j], 4]) / VC[j];
        c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]);
        t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 6]) - t_1;
        
      }else{
        if(solver == 1){
          
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
                           
        }else{
          
          x_ipred[subj_start[j]:subj_end[j],] =
            pmx_solve_bdf(depot_1cmt_linear_plus_mm_ode,
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
                         
        }
        
      }
          
      ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
    
    }

    for(i in 1:n_total){
      if(ipred[i] == 0){
         dv[i] = 0;
      }else{
        dv[i] = lognormal_rng(log(ipred[i]), sigma);
      }
    }
  }
}

