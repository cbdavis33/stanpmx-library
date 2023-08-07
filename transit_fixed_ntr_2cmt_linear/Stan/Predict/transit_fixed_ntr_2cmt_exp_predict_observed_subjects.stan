// Two-compartment PK Model with Transit Compartment absorption
// IIV on CL, VC, Q, VP, KA, MTT (full covariance matrix)
// n_transit is fixed to a positive integer that can be defined as data
// The 0th transit compartment is where the dosing happens (cmt = 1). This could
//   also be called the Depot
// General ODE solution using Torsten to get out individual estimates of 
//   AUC, Cmax, Tmax, ...
// exponential error - DV = IPRED*exp(eps)

functions{
  
  vector transit_fixed_2cmt_ode(real t, vector y, array[] real params,
                                array[] real x_r, array[] int x_i){

    real cl = params[1];
    real vc = params[2];
    real q = params[3];
    real vp = params[4];
    real ka = params[5];
    real mtt = params[6];

    real t_1 = x_r[1];
    real t_2 = x_r[2];

    int n_transit = x_i[1];

    real ktr = (n_transit + 1)/mtt;

    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;

    real slope = ka*y[n_transit + 2] - (ke + k_cp)*y[n_transit + 3] +
                          k_pc*y[n_transit + 4];
    real x = slope > 0 && y[n_transit + 3]/vc > y[n_transit + 7] ? slope/vc : 0;
    real z = t <= t_1 || (slope > 0 && t >= t_1 && t <= t_2) ? 1 : 0;

    vector[n_transit + 8] dydt;

    dydt[1] = -ktr*y[1];                               // 0th transit compartment (depot)
    for(i in 2:(n_transit + 1)){
      dydt[i] = ktr*(y[i - 1] - y[i]);
    }
    dydt[n_transit + 2] = ktr*y[n_transit + 1] - ka*y[n_transit + 2];
    dydt[n_transit + 3] = slope;
    dydt[n_transit + 4] = k_cp*y[n_transit + 3] - k_pc*y[n_transit + 4];
    dydt[n_transit + 5] = y[n_transit + 3];                              // AUC
    dydt[n_transit + 6] = t >= t_1 && t <= t_2 ? y[n_transit + 3] : 0;   // AUC_t_1-t_2
    dydt[n_transit + 7] = x;                                             // C_max
    dydt[n_transit + 8] = z;                                             // t_max

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
  
  int<lower = 1> n_transit;
  
  real<lower = 0> t_1;   // Time at which to start SS calculations (AUC_ss, C_max_ss, ...)
  real<lower = t_1> t_2; // Time at which to end SS calculations (AUC_ss, C_max_ss, ...)
 
}
transformed data{ 
  
  int n_random = 6;           // Number of random effects
  int n_cmt = n_transit + 8;  // Depot, tr_1, ..., tr_n, absorption, central, peripheral
                              //   AUC, AUC_t1-t2, Cmax, Tmax
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt);
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  array[1, 2] real x_r = {{t_1, t_2}};
  
  array[1, 1] int x_i = {{n_transit}};
  
}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVQ;       
  real<lower = 0> TVVP;
  real<lower = 0.5*(TVCL/TVVC + TVQ/TVVC + TVQ/TVVP +
    sqrt((TVCL/TVVC + TVQ/TVVC + TVQ/TVVP)^2 - 4*TVCL/TVVC*TVQ/TVVP))> TVKA;
  real<lower = 0> TVMTT;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{
  
  vector[n_time_new] ipred;   // ipred for the observed individuals at the new timepoints
  vector[n_time_new] pred;    // pred for the observed individuals at the new timepoints
  vector[n_time_new] dv;      // dv for the observed individuals at the new timepoints
  vector[n_time_new] auc;     // auc for the observed individuals from time 0 to the new timepoint
  vector[n_subjects] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[n_subjects] c_max;   // Cmax between t1 and t2 (c_max_ss)
  vector[n_subjects] t_max;   // Tmax between t1 and t2, then subtract off t1
  vector[n_subjects] t_half_alpha;  // alpha half-life
  vector[n_subjects] t_half_terminal;  // terminal half-life
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;
  vector[n_subjects] KA;
  vector[n_subjects] MTT;
  vector[n_subjects] KTR;
 
  {
    
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP,
                                                         TVKA, TVMTT});

    matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
    matrix[n_random, n_random] Omega = quad_form_diag(R, omega);

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));

    real TVKTR = (n_transit + 1)/TVMTT;

    matrix[n_time_new, n_cmt] x_ipred;
    matrix[n_time_new, n_cmt] x_pred;
    
    vector[n_subjects] alpha;
    vector[n_subjects] beta;
    
    CL = col(theta, 1);
    VC = col(theta, 2);
    Q = col(theta, 3);
    VP = col(theta, 4);
    KA = col(theta, 5);
    MTT = col(theta, 6);
    KTR = (n_transit + 1) ./ MTT;
    
    alpha = 0.5*(CL./VC + Q./VC + Q./VP + 
                 sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
    beta = 0.5*(CL./VC + Q./VC + Q./VP - 
                 sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
    
    for(j in 1:n_subjects){
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(transit_fixed_2cmt_ode,
                       n_cmt,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {CL[j], VC[j], Q[j], VP[j], KA[j], MTT[j]},
                       bioav, tlag, x_r, x_i)';
                       
      ipred[subj_start[j]:subj_end[j]] =
        x_ipred[subj_start[j]:subj_end[j], (n_transit + 3)] ./ VC[j];
        
      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(transit_fixed_2cmt_ode,
                       n_cmt,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {TVCL, TVVC, TVQ, TVVP, TVKA, TVMTT},
                       bioav, tlag, x_r, x_i)';

      pred[subj_start[j]:subj_end[j]] = 
        x_pred[subj_start[j]:subj_end[j], (n_transit + 3)] ./ TVVC;
 
      auc[subj_start[j]:subj_end[j]] = 
                    x_ipred[subj_start[j]:subj_end[j], n_transit + 5] ./ VC[j];
      
      auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], n_transit + 6]) / VC[j];
      c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], n_transit + 7]);
      t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], n_transit + 8]) - t_1;
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

