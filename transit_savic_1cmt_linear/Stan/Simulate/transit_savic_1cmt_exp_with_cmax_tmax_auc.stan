// One-compartment PK Model with Transit Compartment absorption
// IIV on CL, VC, KA, NTR, MTT (full covariance matrix)
// n_transit is a positive real number, not fixed, and not necessarily an integer
// General ODE solution using pure Stan code
// exponential error - DV = IPRED*exp(eps)
// The 0th transit compartment is where the dosing happens (cmt = 1). This could
//   also be called the Depot
// Drug is absorbed from the depot through the transit compartments into the 
//   absorption compartment and into the central compartment
//   (depot -> tr1 -> ... -> trn -> absorption -> central)
// Output includes individual Cmax over the whole time period, Tmax between t1 
//   and t2, AUC since 0 for every timepoint, and AUC between t1 and t2 (like a 
//   dosing interval)

functions{
  
  vector transit_1cmt_ode(real t, 
                          vector y,
                          real cl, real vc, real ka, 
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
    
    vector[6] dydt;
    real slope = ka*y[1] -ke*y[2];
    real x = (slope > 0 && y[2]/vc > y[5]) ? slope/vc : 0;
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
  real<lower = 0> TVKA;
  real<lower = 0> TVNTR;
  real<lower = 0> TVMTT;
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_ka;
  real<lower = 0> omega_ntr;
  real<lower = 0> omega_mtt;
  
  corr_matrix[5] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_cl_vc, cor_cl_ka, ... and then construct the 
                     // correlation matrix in transformed data, but it's easy
                     // enough to do in R
  
  real<lower = 0> sigma;
  
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
  
  int n_random = 5; // Number of random effects
  int n_cmt = 6;   // Absorption, central, AUC, AUC_t1-t2, Cmax, Tmax

  vector[n_random] omega = [omega_cl, omega_vc, omega_ka, 
                            omega_ntr, omega_mtt]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);
  
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
  vector[n_subjects] t_half;      // half-life
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  vector[n_subjects] NTR;
  vector[n_subjects] MTT;
  vector[n_subjects] KE;
  
  {
    
    vector[n_random] typical_values = to_vector({TVCL, TVVC, TVKA, 
                                                 TVNTR, TVMTT});
    
    matrix[n_random, n_subjects] eta;   
    matrix[n_subjects, n_random] theta; 
  
    array[n_total] vector[n_cmt] x_ipred;
    
    for(i in 1:n_subjects){
      eta[, i] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                           diag_pre_multiply(omega, L));
    }
    theta = (rep_matrix(typical_values, n_subjects) .* exp(eta))';

    CL = col(theta, 1);
    VC = col(theta, 2);
    KA = col(theta, 3);
    NTR = col(theta, 4);
    MTT = col(theta, 5);
    KE = CL ./ VC;
  
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
          ode_rk45(transit_1cmt_ode, y0, t0, time[(subj_start[j] + 1):subj_end[j]], 
                   CL[j], VC[j], KA[j], NTR[j], MTT[j], bioav[j],
                   n_dose_subj, dosetime_subj, doseamt_subj,
                   t_1, t_2);
      }else if(solver == 2){
        x_ipred[(subj_start[j] + 1):subj_end[j],] = 
          ode_bdf(transit_1cmt_ode, y0, t0, time[(subj_start[j] + 1):subj_end[j]], 
                  CL[j], VC[j], KA[j], NTR[j], MTT[j], bioav[j],
                  n_dose_subj, dosetime_subj, doseamt_subj,
                   t_1, t_2);
      }else if(solver == 3){
        x_ipred[(subj_start[j] + 1):subj_end[j],] = 
          ode_adams(transit_1cmt_ode, y0, t0, time[(subj_start[j] + 1):subj_end[j]], 
                    CL[j], VC[j], KA[j], NTR[j], MTT[j], bioav[j],
                    n_dose_subj, dosetime_subj, doseamt_subj,
                   t_1, t_2);
      }else{
        x_ipred[(subj_start[j] + 1):subj_end[j],] = 
          ode_ckrk(transit_1cmt_ode, y0, t0, time[(subj_start[j] + 1):subj_end[j]], 
                   CL[j], VC[j], KA[j], NTR[j], MTT[j], bioav[j],
                   n_dose_subj, dosetime_subj, doseamt_subj,
                   t_1, t_2);
      }
      
      for(k in subj_start[j]:subj_end[j]){
        ipred[k] = x_ipred[k, 2] / VC[j];
        auc[k] = x_ipred[k, 3] / VC[j];
      }
      
      auc_t1_t2[j] = max(x_ipred[subj_start[j]:subj_end[j], 4]) / VC[j];
      c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]);
      t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 6]) - t_1;
      t_half[j] = log(2)/KE[j];
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
