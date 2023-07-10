// One-compartment PK Model with Transit Compartment absorption
// IIV on CL, VC, KA, NTR, MTT (full covariance matrix)
// n_transit is a positive real number, not fixed, and not necessarily an integer
// General ODE solution using pure Stan code to get out individual estimates of 
//   AUC, Cmax, Tmax, ...
// proportional error - DV = IPRED*(1 + eps_p)
// Predictions are generated from a normal that is truncated below at 0

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  vector transit_1cmt_ode(real t, 
                          vector y,
                          real ke, real vc, real ka, real n_tr, real f, 
                          real ktr, 
                          int n_dose_subj,            // # of doses for this subject
                          array[] real dosetime_subj, // dosetimes for this subject
                          array[] real doseamt_subj,  // dose amounts for this subject
                          real t_1, real t_2){ 
    
    // real ktr = (n_tr + 1)/mtt;
    real k_inpt = f*pow(ktr, n_tr + 1)/exp(lgamma(n_tr + 1));
    
    real inpt = 0;
    array[n_dose_subj] real ipt = rep_array(0.0, n_dose_subj);
    
    vector[6] dydt;
    real slope = ka*y[1] - ke*y[2];
    real x = (slope > 0 && y[2]/vc > y[5]) ? slope/vc : 0;
    real z = t <= t_1 || (slope > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    for(i in 1:n_dose_subj){
      if(t >= dosetime_subj[i]){
        real delta_t = t - dosetime_subj[i];
        ipt[i] = doseamt_subj[i]*pow(delta_t, n_tr)*exp(-ktr*delta_t);
      }
    }
    
    inpt = sum(ipt);
    
    dydt[1] = k_inpt*inpt - ka*y[1];
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
  
  int n_dose;
  array[n_dose] real dosetime;
  array[n_dose] real doseamt;
  
  array[n_subjects_new] int subj_start_dose;
  array[n_subjects_new] int subj_end_dose;
  
  real<lower = 0> t_1;   // Time at which to start SS calculations (AUC_ss, C_max_ss, ...)
  real<lower = t_1> t_2; // Time at which to end SS calculations (AUC_ss, C_max_ss, ...)
  
}
transformed data{ 
  
  int n_random = 5; // Number of random effects
  int n_cmt = 6;   // Absorption, central, AUC, AUC_t1-t2, Cmax, Tmax
  
  array[n_subjects] real bioav = rep_array(1.0, n_subjects);
  
  vector[n_cmt] y0 = rep_vector(0.0, n_cmt);

}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = TVCL/TVVC> TVKA;
  real<lower = 0> TVNTR; 
  real<lower = 0> TVMTT; 
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{

  vector[n_time_new] ipred;   // ipred for the observed individuals at the new timepoints
  vector[n_time_new] pred;    // pred for the observed individuals at the new timepoints
  vector[n_time_new] dv;      // dv for the observed individuals at the new timepoints
  vector[n_time_new] auc;     // auc for the observed individuals from time 0 to the new timepoint
  vector[n_subjects_new] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[n_subjects_new] c_max;   // Cmax between t1 and t2 (c_max_ss)
  vector[n_subjects_new] t_max;   // Tmax between t1 and t2, then subtract off t1
  vector[n_subjects_new] t_half;  // half-life
  
  vector[n_subjects_new] CL;
  vector[n_subjects_new] VC;
  vector[n_subjects_new] KA;
  vector[n_subjects_new] NTR;
  vector[n_subjects_new] MTT;
  vector[n_subjects_new] KE;
  vector[n_subjects_new] KTR;

  {
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA, 
                                                         TVNTR, TVMTT});
    
    matrix[n_subjects_new, n_random] eta_new;
    matrix[n_subjects_new, n_random] theta_new;
    
    array[n_time_new] vector[n_cmt] x_pred;
    array[n_time_new] vector[n_cmt] x_ipred;
    
    real TVKE = TVCL/TVVC;
    real TVKTR = (TVNTR + 1)/TVMTT;
    
    for(i in 1:n_subjects_new){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects_new) .* exp(eta_new));
    
    CL = col(theta_new, 1);
    VC = col(theta_new, 2);
    KA = col(theta_new, 3);
    NTR = col(theta_new, 4);
    MTT = col(theta_new, 5);
    KE = CL ./ VC;
    KTR = (NTR + 1) ./ MTT;

    for(j in 1:n_subjects_new){
      
      real t0 = time[subj_start[j]];
      int n_dose_subj = subj_end_dose[j] - subj_start_dose[j] + 1;
      
      array[n_dose_subj] real dosetime_subj = 
        dosetime[subj_start_dose[j]:subj_end_dose[j]];
      array[n_dose_subj] real doseamt_subj = 
        doseamt[subj_start_dose[j]:subj_end_dose[j]];
      
      x_ipred[subj_start[j],] = y0;
      x_pred[subj_start[j],] = y0;
      
      x_ipred[(subj_start[j] + 1):subj_end[j],] =
          ode_ckrk(transit_1cmt_ode, y0, t0,
                   time[(subj_start[j] + 1):subj_end[j]],
                   KE[j], VC[j], KA[j], NTR[j], bioav[j], KTR[j],
                   n_dose_subj, dosetime_subj, doseamt_subj,
                   t_1, t_2);
                      
      x_pred[(subj_start[j] + 1):subj_end[j],] = 
          ode_ckrk(transit_1cmt_ode, y0, t0, 
                   time[(subj_start[j] + 1):subj_end[j]], 
                   TVKE, TVVC, TVKA, TVNTR, bioav[j], TVKTR,
                   n_dose_subj, dosetime_subj, doseamt_subj,
                   t_1, t_2);

      for(k in subj_start[j]:subj_end[j]){
        ipred[k] = fmax(0.000000001, x_ipred[k, 2] / VC[j]);
        pred[k] = fmax(0.000000001, x_pred[k, 2] / TVVC);
        auc[k] = x_ipred[k, 3] / VC[j];
      }
 
      auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], 4]) / VC[j];
      c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]);
      t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 6]) - t_1;
      t_half[j] = log(2)/KE[j];
    }

    for(i in 1:n_time_new){
      if(ipred[i] == 0){
        dv[i] = 0;
      }else{
        real ipred_tmp = ipred[i];
        real sigma_tmp = ipred_tmp*sigma_p;
        dv[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
      }
    }
  
  }

}


