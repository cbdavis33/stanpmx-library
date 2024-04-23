// First Order Absorption (oral/subcutaneous)
// Two-compartment PK Model, Friberg-Karlsson neutropenia PD model
// IIV on CL, VC, Q, VP, and Ka (full covariance matrix)
// IIV on MTT, CIRC0, GAMMA, ALPHA (full covariance matrix) for PD
// proportional plus additive error for PK - DV = IPRED*(1 + eps_p) + eps_a
// proportional plus additive error for PD - DV = IPRED*(1 + eps_p_pd) + eps_a_pd
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
  
  vector depot_2cmt_friberg_ode(real t, vector y, array[] real params, 
                                array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real q = params[3];
    real vp = params[4];
    real ka = params[5];
    real mtt = params[6];
    real circ_0 = params[7];
    real gamma = params[8];
    real alpha = params[9];
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;
    real k_tr = 4/mtt; // k_tr = (n_tr + 1)/mtt    
    
    real conc = y[2]/vc;
    
    real e_drug = fmin(1.0, alpha*conc); // Maybe reparameterize this so no more fmin?
    real prol = y[4] + circ_0;
    real transit_1 = y[5] + circ_0; 
    real transit_2 = y[6] + circ_0;
    real transit_3 = y[7] + circ_0;
    real circ = y[8] + circ_0; // fmax(machine_precision(), y[8] + circ_0)
    
    real slope_pk = ka*y[1] - (ke + k_cp)*y[2] + k_pc*y[3];
    real x_pk = slope_pk > 0 && conc > y[11] ? slope_pk/vc : 0;
    real z_pk = t <= t_1 || (slope_pk > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    real slope_pd = k_tr*(transit_3 - circ);
    real x_pd = slope_pd < 0 && y[8] < y[13] ? slope_pd : 0;
    
    vector[13] dydt;
    
    dydt[1] = -ka*y[1];                               // depot
    dydt[2] = ka*y[1] - (ke + k_cp)*y[2] + k_pc*y[3]; // central
    dydt[3] = k_cp*y[2] - k_pc*y[3];                  // peripheral
    dydt[4] = k_tr*prol*((1 - e_drug)*(circ_0/circ)^gamma - 1);  // proliferative cells
    dydt[5] = k_tr*(prol - transit_1);                // transit 1
    dydt[6] = k_tr*(transit_1 - transit_2);           // transit 2
    dydt[7] = k_tr*(transit_2 - transit_3);           // transit 3
    dydt[8] = slope_pd;                               // circulating blood cells
    dydt[9] = y[2];                                   // AUC
    dydt[10] = t >= t_1 && t <= t_2 ? y[2] : 0;       // AUC_t_1-t_2
    dydt[11] = x_pk;                                  // C_max 
    dydt[12] = z_pk;                                  // t_max for PK
    dydt[13] = x_pd;                                  // R_min
    
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
  
  int n_random = 5;
  int n_random_pd = 4;
  
  int n_cmt = 3;       // number of ODEs in PK model (depot, central, peripheral)
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
  real<lower = 0> TVQ;       
  real<lower = 0> TVVP;
  real<lower = 0.5*(TVCL/TVVC + TVQ/TVVC + TVQ/TVVP +
    sqrt((TVCL/TVVC + TVQ/TVVC + TVQ/TVVP)^2 - 4*TVCL/TVVC*TVQ/TVVP))> TVKA;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  vector<lower = 0>[2] sigma;
  cholesky_factor_corr[2] L_Sigma;
  
  matrix[n_random, n_subjects] Z;
  
  real<lower = 0> TVMTT;       
  real<lower = 0> TVCIRC0; 
  real<lower = 0> TVGAMMA;
  real<lower = 0> TVALPHA;
  
  vector<lower = 0>[n_random_pd] omega_pd;
  cholesky_factor_corr[n_random_pd] L_pd;
  
  vector<lower = 0>[2] sigma_pd;
  cholesky_factor_corr[2] L_Sigma_pd;
  
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
  vector[n_subjects_new] t_half_alpha;    // alpha half-life
  vector[n_subjects_new] t_half_terminal; // terminal half-life
  vector[n_subjects_new] r_min;       // Minimum response
  
  vector[n_subjects_new] CL;
  vector[n_subjects_new] VC;
  vector[n_subjects_new] Q;
  vector[n_subjects_new] VP;
  vector[n_subjects_new] KA;
  
  vector[n_subjects_new] MTT;
  vector[n_subjects_new] CIRC0;
  vector[n_subjects_new] GAMMA;
  vector[n_subjects_new] ALPHA;

  {
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP, 
                                                         TVKA});
    
    matrix[n_subjects_new, n_random] eta_new;
    matrix[n_subjects_new, n_random] theta_new;
    
    row_vector[n_random_pd] typical_values_pd = to_row_vector({TVMTT, TVCIRC0, 
                                                               TVGAMMA, TVALPHA});
    
    matrix[n_subjects_new, n_random_pd] eta_new_pd;   
    matrix[n_subjects_new, n_random_pd] theta_new_pd; 
    
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
    
    real sigma_sq_p = Sigma[1, 1];
    real sigma_sq_a = Sigma[2, 2];
    real sigma_p_a = Sigma[1, 2];
    
    matrix[2, 2] R_Sigma_pd = multiply_lower_tri_self_transpose(L_Sigma_pd);
    matrix[2, 2] Sigma_pd = quad_form_diag(R_Sigma_pd, sigma_pd);
    
    real sigma_sq_p_pd = Sigma_pd[1, 1];
    real sigma_sq_a_pd = Sigma_pd[2, 2];
    real sigma_p_a_pd = Sigma_pd[1, 2];
    
    vector[n_subjects_new] alpha;
    vector[n_subjects_new] beta;
    
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
    Q = col(theta_new, 3);
    VP = col(theta_new, 4);
    KA = col(theta_new, 5);
    
    MTT = col(theta_new_pd, 1);
    CIRC0 = col(theta_new_pd, 2);
    GAMMA = col(theta_new_pd, 3);
    ALPHA = col(theta_new_pd, 4);
    
    alpha = 0.5*(CL./VC + Q./VC + Q./VP + 
                 sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
    beta = 0.5*(CL./VC + Q./VC + Q./VP - 
                 sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));

    for(j in 1:n_subjects_new){
      
      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_2cmt_friberg_ode,
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
                        MTT[j], CIRC0[j], GAMMA[j], ALPHA[j]},
                       bioav, tlag, x_r)';

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_2cmt_friberg_ode,
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
                        TVMTT, TVCIRC0, TVGAMMA, TVALPHA},
                       bioav, tlag, x_r)';

      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){
          ipred[k] = x_ipred[k, 2] / VC[j];
          pred[k] = x_pred[k, 2] / TVVC;
        }else if(cmt[k] == 4){
          ipred[k] = x_ipred[k, 8] + CIRC0[j];
          pred[k] = x_pred[k, 8] + TVCIRC0;
        }
      }
      
      auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 9] ./ VC[j];  
      auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], 10]) / VC[j];
      c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 11]);
      t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 12]) - t_1;
      t_half_alpha[j] = log(2)/alpha[j];
      t_half_terminal[j] = log(2)/beta[j];
      
      r_min[j] = min(x_ipred[subj_start[j]:subj_end[j], 13]) + CIRC0[j];
    
    }

    for(i in 1:n_time_new){
      if(ipred[i] == 0){
        dv[i] = 0;
      }else if(cmt[i] == 2 || cmt[i] == 4){
        real ipred_tmp = ipred[i];
        real sigma_tmp = cmt[i] == 2 ? 
        sqrt(square(ipred_tmp) * sigma_sq_p + sigma_sq_a + 2*ipred_tmp*sigma_p_a) : 
        sqrt(square(ipred_tmp) * sigma_sq_p_pd + sigma_sq_a_pd + 
                          2*ipred_tmp*sigma_p_a_pd);
    
        dv[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
      }
    }
  }
}
