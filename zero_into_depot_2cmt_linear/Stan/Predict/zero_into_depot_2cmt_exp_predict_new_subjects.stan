// First Order Absorption (oral/subcutaneous)
// Two-compartment PK Model
// Zero-order distributive delay (like an infusion into the gut) to bring about 
//   a delayed absorption
// IIV on CL, VC, Q, VP, KA, DUR (full covariance matrix)
// exponential error - DV = IPRED*exp(eps)
// General ODE solution using Torsten to get out individual estimates of AUC, 
//   Cmax, Tmax, ...

functions{
  
  vector depot_2cmt_ode(real t, vector y, array[] real params, 
                        array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real q = params[3];
    real vp = params[4];
    real ka = params[5];
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;
    
    real slope = ka*y[1] - (ke + k_cp)*y[2] + k_pc*y[3];
    real x = slope > 0 && y[2]/vc > y[6] ? slope/vc : 0;
    real z = t <= t_1 || (slope > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    vector[7] dydt;

    dydt[1] = -ka*y[1];                            // depot
    dydt[2] = slope;                               // central
    dydt[3] = k_cp*y[2] - k_pc*y[3];               // peripheral
    dydt[4] = y[2];                                // AUC
    dydt[5] = t >= t_1 && t <= t_2 ? y[2] : 0;     // AUC_t_1-t_2
    dydt[6] = x;                                   // C_max
    dydt[7] = z;                                   // t_max
    
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
  array[n_time_new] real ii;
  array[n_time_new] int addl;
  array[n_time_new] int ss;
  array[n_subjects_new] int subj_start;
  array[n_subjects_new] int subj_end;
  
  real<lower = 0> t_1;   // Time at which to start SS calculations (AUC_ss, C_max_ss, ...)
  real<lower = t_1> t_2; // Time at which to end SS calculations (AUC_ss, C_max_ss, ...)
  
}
transformed data{ 
  
  int n_random = 6; // Number of random effects
  int n_cmt = 7;    // Number of compartments (depot, central, peripheral, AUC, AUC_ss, Cmax_ss, Tmax_ss)
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt);
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  array[1, 2] real x_r = {{t_1, t_2}};

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
  
  real<lower = 0> sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{

  vector[n_time_new] ipred;       // ipred for the observed individuals at the new timepoints
  vector[n_time_new] pred;        // pred for the observed individuals at the new timepoints
  vector[n_time_new] dv;          // dv for the observed individuals at the new timepoints
  vector[n_time_new] auc;         // auc for the observed individuals from time 0 to the new timepoint
  vector[n_subjects_new] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[n_subjects_new] c_max;   // Cmax between t1 and t2 (c_max_ss)
  vector[n_subjects_new] t_max;   // Tmax between t1 and t2, then subtract off t1
  vector[n_subjects_new] t_half_alpha;    // alpha half-life
  vector[n_subjects_new] t_half_terminal; // terminal half-life
  
  vector[n_subjects_new] CL;
  vector[n_subjects_new] VC;
  vector[n_subjects_new] Q;
  vector[n_subjects_new] VP;
  vector[n_subjects_new] KA;
  vector[n_subjects_new] DUR;
  vector[n_subjects_new] KE;

  {
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP,
                                                        TVKA, TVDUR});
    
    matrix[n_subjects_new, n_random] eta_new;
    matrix[n_subjects_new, n_random] theta_new;
    matrix[n_time_new, n_cmt] x_ipred;
    matrix[n_time_new, 3] x_pred;
    
    vector[n_subjects_new] alpha;
    vector[n_subjects_new] beta;

    for(i in 1:n_subjects_new){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects_new) .* exp(eta_new));
    
    CL = col(theta_new, 1);
    VC = col(theta_new, 2);
    Q = col(theta_new, 3);
    VP = col(theta_new, 4);
    KA = col(theta_new, 5);
    DUR = col(theta_new, 6);
    
    alpha = 0.5*(CL./VC + Q./VC + Q./VP + 
                 sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
    beta = 0.5*(CL./VC + Q./VC + Q./VP - 
                 sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));

    for(j in 1:n_subjects_new){

      array[subj_end[j] - subj_start[j] + 1] real rate = 
              to_array_1d(to_vector(amt[subj_start[j]:subj_end[j]]) ./ DUR[j]);
              
      array[subj_end[j] - subj_start[j] + 1] real rate_tv = 
              to_array_1d(to_vector(amt[subj_start[j]:subj_end[j]]) ./ TVDUR);
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_2cmt_ode,
                       n_cmt,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate,
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {CL[j], VC[j], Q[j], VP[j], KA[j]}, 
                       bioav, tlag, x_r)';
                      
      ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate_tv,
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {TVCL, TVQ, TVVC, TVVP, TVKA})';

      pred[subj_start[j]:subj_end[j]] = 
        x_pred[subj_start[j]:subj_end[j], 2] ./ TVVC;
        
      auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 4] ./ VC[j];    
      auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]) / VC[j];
      c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 6]);
      t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 7]) - t_1;
      t_half_alpha[j] = log(2)/alpha[j];
      t_half_terminal[j] = log(2)/beta[j];
    
    }

    for(i in 1:n_time_new){
      if(ipred[i] == 0){
        dv[i] = 0;
      }else{
        dv[i] = lognormal_rng(log(ipred[i]), sigma);
      }
    }
  
  }

}
