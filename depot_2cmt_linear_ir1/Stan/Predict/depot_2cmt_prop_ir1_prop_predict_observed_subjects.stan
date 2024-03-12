// First Order Absorption (oral/subcutaneous)
// Two-compartment PK Model, IR1 PD model
// IIV on CL, VC, Q, VP, and Ka (full covariance matrix)
// IIV on KIN, KOUT, IC50 (full covariance matrix) for PD. IMAX and HILL are 
//   fixed to be 1
// proportional error for PK - DV = IPRED*(1 + eps_p)
// proportional error for PD - DV = IPRED*(1 + eps_p_pd)
// General ODE solution using Torsten to get out individual estimates of AUC, 
//   Cmax, Tmax, ...
// Predictions are generated from a normal that is truncated below at 0

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  vector depot_2cmt_ir1_ode(real t, vector y, array[] real params, 
                            array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real q = params[3];
    real vp = params[4];
    real ka = params[5];
    real kin = params[6];
    real kout = params[7];
    real ic50 = params[8];
    real imax = params[9];   // It's fixed to 1 in this particular model
    real hill = params[10];  // It's fixed to 1 in this particular model
    real r_0 = params[11];
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;
    
    real conc = y[2]/vc;
    
    real inh = (imax*pow(conc, hill))/(pow(ic50, hill) + pow(conc, hill));
    real response = y[4] + r_0;
    
    real slope_pk = ka*y[1] - (ke + k_cp)*y[2] + k_pc*y[3];
    real x_pk = slope_pk > 0 && conc > y[7] ? slope_pk/vc : 0;
    real z_pk = t <= t_1 || (slope_pk > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    real slope_pd = kin*(1 - inh) - kout*response;
    real x_pd = slope_pd < 0 && y[4] < y[9] ? slope_pd : 0;
    real z_pd = t <= t_1 || (slope_pd < 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    vector[10] dydt;

    dydt[1] = -ka*y[1];
    dydt[2] = slope_pk;
    dydt[3] = k_cp*y[2] - k_pc*y[3];
    dydt[4] = slope_pd;
    dydt[5] = y[2];                                // AUC
    dydt[6] = t >= t_1 && t <= t_2 ? y[2] : 0;     // AUC_t_1-t_2
    dydt[7] = x_pk;                                // C_max 
    dydt[8] = z_pk;                                // t_max for PK
    dydt[9] = x_pd;                                // R_min
    dydt[10] = z_pd;                               // t_min for PD
    
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
  
  int n_random = 5;
  int n_random_pd = 3;
  
  int n_cmt = 3;       // number of ODEs in PK model (depot, central, peripheral)
  int n_cmt_pd = 1;    // number of ODEs in PD system
  int n_cmt_extra = 6; // number of ODEs for AUC, Tmax, Rmin, ...
  
  array[n_cmt + n_cmt_pd + n_cmt_extra] real bioav = 
                                 rep_array(1.0, n_cmt + n_cmt_pd + n_cmt_extra); 
  array[n_cmt + n_cmt_pd + n_cmt_extra] real tlag = 
                                 rep_array(0.0, n_cmt + n_cmt_pd + n_cmt_extra);
  
  real TVIMAX = 1.0;
  real TVHILL = 1.0;
  vector<lower = 0>[n_subjects] IMAX = rep_vector(TVIMAX, n_subjects); // IMAX and HILL are both fixed to 1.0 in this model,
  vector<lower = 0>[n_subjects] HILL = rep_vector(TVHILL, n_subjects); // but it could be data or a parameter in another model.
             
  
  array[1, 2] real x_r = {{t_1, t_2}};

}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVQ;       
  real<lower = 0> TVVP;
  real<lower = 0.5*(TVCL/TVVC + TVQ/TVVC + TVQ/TVVP +
    sqrt((TVCL/TVVC + TVQ/TVVC + TVQ/TVVP)^2 - 4*TVCL/TVVC*TVQ/TVVP))> TVKA;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
  real<lower = 0> TVKIN;       
  real<lower = 0> TVKOUT; 
  real<lower = 0> TVIC50;
  
  vector<lower = 0>[n_random_pd] omega_pd;
  cholesky_factor_corr[n_random_pd] L_pd;
  
  real<lower = 0> sigma_p_pd;
  
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
  vector[n_subjects] t_half_alpha;    // alpha half-life
  vector[n_subjects] t_half_terminal; // terminal half-life
  vector[n_subjects] r_min;       // Minimum response
  vector[n_subjects] t_min;       // Tmin for PD
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;
  vector[n_subjects] KA;
  
  vector[n_subjects] KIN;
  vector[n_subjects] KOUT;
  vector[n_subjects] IC50;
 
  {
    
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP, 
                                                         TVKA});

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
                          
    row_vector[n_random_pd] typical_values_pd = to_row_vector({TVKIN, TVKOUT, 
                                                               TVIC50});
    
    matrix[n_subjects, n_random_pd] eta_pd = 
                                      diag_pre_multiply(omega_pd, L_pd * Z_pd)';  
    matrix[n_subjects, n_random_pd] theta_pd =
                     (rep_matrix(typical_values_pd, n_subjects) .* exp(eta_pd));
    
    matrix[n_time_new, n_cmt + n_cmt_pd + n_cmt_extra] x_pred;
    matrix[n_time_new, n_cmt + n_cmt_pd + n_cmt_extra] x_ipred;
    
    vector[n_subjects] r_0;
    vector[n_subjects] alpha;
    vector[n_subjects] beta;
    
    real TVR0 = TVKIN/TVKOUT;
    
    CL = col(theta, 1);
    VC = col(theta, 2);
    Q = col(theta, 3);
    VP = col(theta, 4);
    KA = col(theta, 5);
    
    KIN = col(theta_pd, 1);
    KOUT = col(theta_pd, 2);
    IC50 = col(theta_pd, 3);
    
    r_0 = KIN ./ KOUT;
    alpha = 0.5*(CL./VC + Q./VC + Q./VP + 
                 sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
    beta = 0.5*(CL./VC + Q./VC + Q./VP - 
                 sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));

    for(j in 1:n_subjects){
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_2cmt_ir1_ode,
                       n_cmt + n_cmt_pd + n_cmt_extra,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {CL[j], VC[j], Q[j], VP[j], KA[j], 
                        KIN[j], KOUT[j], IC50[j], TVIMAX, TVHILL, r_0[j]},
                       bioav, tlag, x_r)';

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_2cmt_ir1_ode,
                       n_cmt + n_cmt_pd + n_cmt_extra,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {TVCL, TVVC, TVQ, TVVP, TVKA, 
                        TVKIN, TVKOUT, TVIC50, TVIMAX, TVHILL, TVR0},
                       bioav, tlag, x_r)';

      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){
          ipred[k] = x_ipred[k, 2] / VC[j];
          pred[k] = x_pred[k, 2] / TVVC;
        }else if(cmt[k] == 4){
          ipred[k] = x_ipred[k, 4] + r_0[j];
          pred[k] = x_pred[k, 4] + TVR0;
        }
      }
      
      auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 5] ./ VC[j];  
      auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], 6]) / VC[j];
      c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 7]);
      t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 8]) - t_1;
      t_half_alpha[j] = log(2)/alpha[j];
      t_half_terminal[j] = log(2)/beta[j];
      
      r_min[j] = min(x_ipred[subj_start[j]:subj_end[j], 9]) + r_0[j];
      t_min[j] = max(x_ipred[subj_start[j]:subj_end[j], 10]) - t_1;
    }
    
    for(i in 1:n_time_new){
      if(ipred[i] == 0){
        dv[i] = 0;
      }else if(cmt[i] == 2 || cmt[i] == 4){
        real ipred_tmp = ipred[i];
        real sigma_tmp = cmt[i] == 2 ? ipred_tmp*sigma_p : ipred_tmp*sigma_p_pd;
    
        dv[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
      }
    }
  }
}
