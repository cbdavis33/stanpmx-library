// Two-compartment PK Model with Transit Compartment absorption
// IIV on CL, VC, Q, VP, KA, NTR, MTT (full covariance matrix)
// n_transit is a positive real number, not fixed, and not necessarily an integer
// General ODE solution using pure Stan code
// proportional plus additive error - DV = IPRED(1 + eps_p) + eps_a
// The 0th transit compartment is where the dosing happens (cmt = 1). This could
//   also be called the Depot
// Drug is absorbed from the depot through the transit compartments into the 
//   absorption compartment and into the central compartment
//   (depot -> tr1 -> ... -> trn -> absorption -> central)
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   generates values from a normal that is truncated below at 0
// Output includes individual Cmax over the whole time period, Tmax between t1 
//   and t2, AUC since 0 for every timepoint, and AUC between t1 and t2 (like a 
//   dosing interval)

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  vector transit_2cmt_ode(real t, 
                          vector y,
                          real cl, real vc, real q, real vp, real ka, 
                          real n_tr, real mtt, real f,
                          int n_dose_subj,             // # of doses for this subject
                          array[] real dosetime_subj,  // dosetimes for this subject
                          array[] real doseamt_subj,   // dose amounts for this subject
                          real t_1, real t_2){ 
    
    real ktr = (n_tr + 1)/mtt;
    // real k_inpt = f*pow(ktr, n_tr + 1)/exp(lgamma(n_tr + 1));
    real log_k_inpt = log(f) + (n_tr + 1)*log(ktr) - lgamma(n_tr + 1);
    real log_inpt;
    array[n_dose_subj] real log_ipt = rep_array(negative_infinity(), n_dose_subj); 
    
    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;
    
    vector[7] dydt;
    real slope = ka*y[1] - (ke + k_cp)*y[2] + k_pc*y[3];
    real x = (slope > 0 && y[2]/vc > y[6]) ? slope/vc : 0;
    real z = t <= t_1 || (slope > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    for(i in 1:n_dose_subj){
      if(t >= dosetime_subj[i]){
        real delta_t = t - dosetime_subj[i];
        // ipt[i] = doseamt_subj[i]*pow(delta_t, n_tr)*exp(-ktr*delta_t);
        log_ipt[i] = log(doseamt_subj[i]) + n_tr*log(delta_t) - ktr*delta_t;
      }
    }
    
    // inpt = sum(ipt);
    // inpt = exp(log_sum_exp(log_ipt));
    log_inpt = log_sum_exp(log_ipt);
    
    dydt[1] = exp(log_k_inpt + log_inpt) - ka*y[1];
    dydt[2] = slope;
    dydt[3] = k_cp*y[2] - k_pc*y[3];
    dydt[4] = y[2];                                // AUC
    dydt[5] = t >= t_1 && t <= t_2 ? y[2] : 0;     // AUC_t_1-t_2
    dydt[6] = x;                                   // C_max 
    dydt[7] = z;                                   // t_max 

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
  real<lower = 0> TVQ;
  real<lower = 0> TVVP;
  real<lower = 0> TVKA;
  real<lower = 0> TVNTR;
  real<lower = 0> TVMTT;
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_q;
  real<lower = 0> omega_vp;
  real<lower = 0> omega_ka;
  real<lower = 0> omega_ntr;
  real<lower = 0> omega_mtt;
  
  corr_matrix[7] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_cl_vc, cor_cl_q, ... and then construct the 
                     // correlation matrix in transformed data, but it's easy
                     // enough to do in R
  
  real<lower = 0> sigma_p;
  real<lower = 0> sigma_a;
  real<lower = -1, upper = 1> cor_p_a;
  
  int n_dose;
  array[n_dose] real dosetime;
  array[n_dose] real doseamt;
  
  array[n_subjects] int subj_start_dose;
  array[n_subjects] int subj_end_dose;
  
  real<lower = 0> t_1;
  real<lower = 0> t_2;
  
  int<lower = 1, upper = 4> solver; // 1 = rk45, 2 = bdf, 3 = adams, 4 = ckrk
  
}
transformed data{
  
  int n_random = 7; // Number of random effects
  int n_cmt = 7;   // Absorption, central, peripheral, AUC, AUC_t1-t2, Cmax, Tmax

  vector[n_random] omega = [omega_cl, omega_vc, omega_q, omega_vp, omega_ka,
                            omega_ntr, omega_mtt]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);

  vector[2] sigma = [sigma_p, sigma_a]';
  matrix[2, 2] R_Sigma = rep_matrix(1, 2, 2);
  R_Sigma[1, 2] = cor_p_a;
  R_Sigma[2, 1] = cor_p_a;
  
  matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
  
  array[n_subjects] real bioav = rep_array(1.0, n_subjects); // Hardcoding, but could be data or a parameter in another situation
  
  vector[n_cmt] y0 = rep_vector(0.0, n_cmt);
  
}
model{
  
}
generated quantities{
  
  vector[n_total] ipred;          // concentration with no residual error
  vector[n_total] dv;             // concentration with residual error
  
  vector[n_total] auc;            // AUC
  vector[n_subjects] auc_t1_t2;   // AUC from t1 up to t2
  vector[n_subjects] c_max;       // Cmax
  vector[n_subjects] t_max;       // Tmax
  vector[n_subjects] t_half_alpha;    // alpha half-life
  vector[n_subjects] t_half_terminal; // terminal half-life
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;
  vector[n_subjects] KA;
  vector[n_subjects] NTR;
  vector[n_subjects] MTT;
  
  {
    
    vector[n_random] typical_values = to_vector({TVCL, TVVC, TVQ, TVVP, 
                                                 TVKA, TVNTR, TVMTT});
    
    matrix[n_random, n_subjects] eta;   
    matrix[n_subjects, n_random] theta; 
  
    array[n_total] vector[n_cmt] x_ipred;
    
    vector[n_subjects] alpha;
    vector[n_subjects] beta;
    
    for(i in 1:n_subjects){
      eta[, i] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                           diag_pre_multiply(omega, L));
    }
    theta = (rep_matrix(typical_values, n_subjects) .* exp(eta))';

    CL = col(theta, 1);
    VC = col(theta, 2);
    Q = col(theta, 3);
    VP = col(theta, 4);
    KA = col(theta, 5);
    NTR = col(theta, 6);
    MTT = col(theta, 7);
    
    alpha = 0.5*(CL./VC + Q./VC + Q./VP + 
                 sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
    beta = 0.5*(CL./VC + Q./VC + Q./VP - 
                 sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
  
    for(j in 1:n_subjects){
      
      real t0 = time[subj_start[j]];
      int n_dose_subj = subj_end_dose[j] - subj_start_dose[j] + 1;
      
      array[n_dose_subj] real dosetime_subj = 
        dosetime[subj_start_dose[j]:subj_end_dose[j]];
      array[n_dose_subj] real doseamt_subj = 
        doseamt[subj_start_dose[j]:subj_end_dose[j]];
      
      x_ipred[subj_start[j],] = y0;
      
      if(solver == 1){
        x_ipred[(subj_start[j] + 1):subj_end[j],] = 
          ode_rk45(transit_2cmt_ode, y0, t0, time[(subj_start[j] + 1):subj_end[j]], 
                   CL[j], VC[j], Q[j], VP[j], KA[j], NTR[j], MTT[j], bioav[j],
                   n_dose_subj, dosetime_subj, doseamt_subj,
                   t_1, t_2);
      }else if(solver == 2){
        x_ipred[(subj_start[j] + 1):subj_end[j],] = 
          ode_bdf(transit_2cmt_ode, y0, t0, time[(subj_start[j] + 1):subj_end[j]], 
                  CL[j], VC[j], Q[j], VP[j], KA[j], NTR[j], MTT[j], bioav[j],
                  n_dose_subj, dosetime_subj, doseamt_subj,
                   t_1, t_2);
      }else if(solver == 3){
        x_ipred[(subj_start[j] + 1):subj_end[j],] = 
          ode_adams(transit_2cmt_ode, y0, t0, time[(subj_start[j] + 1):subj_end[j]], 
                    CL[j], VC[j], Q[j], VP[j], KA[j], NTR[j], MTT[j], bioav[j],
                    n_dose_subj, dosetime_subj, doseamt_subj,
                   t_1, t_2);
      }else{
        x_ipred[(subj_start[j] + 1):subj_end[j],] = 
          ode_ckrk(transit_2cmt_ode, y0, t0, time[(subj_start[j] + 1):subj_end[j]], 
                   CL[j], VC[j], Q[j], VP[j], KA[j], NTR[j], MTT[j], bioav[j],
                   n_dose_subj, dosetime_subj, doseamt_subj,
                   t_1, t_2);
      }
      
      for(k in subj_start[j]:subj_end[j]){
        ipred[k] = x_ipred[k, 2] / VC[j];
        auc[k] = x_ipred[k, 4] / VC[j];
      }
      
      auc_t1_t2[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]) / VC[j];
      c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 6]);
      t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 7]) - t_1;
      t_half_alpha[j] = log(2)/alpha[j];
      t_half_terminal[j] = log(2)/beta[j];
    }
    
    for(i in 1:n_total){
      if(ipred[i] == 0){
         dv[i] = 0;
      }else{
        real ipred_tmp = ipred[i];
        real sigma_tmp = sqrt(square(ipred_tmp) * Sigma[1, 1] + Sigma[2, 2] + 
                              2*ipred_tmp*Sigma[2, 1]);
        dv[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
        
      }
    }
  }
}

