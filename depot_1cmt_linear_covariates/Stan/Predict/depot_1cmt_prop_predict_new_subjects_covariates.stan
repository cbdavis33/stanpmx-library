// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// IIV on CL, VC, and Ka (full covariance matrix)
// proportional error - DV = IPRED*(1 + eps_p)
// User's choice of analytical solution or general ODE solution
// General ODE solution using Torsten will get out individual estimates of AUC, 
//   Cmax, Tmax, ... Analytical will not
// Predictions are generated from a normal that is truncated below at 0
// Covariates: 
//   1) Body Weight on CL and VC - (wt/70)^theta
//   2) Concomitant administration of protein pump inhibitors (CMPPI) 
//      on KA (0/1) - exp(theta*cmppi)
//   3) eGFR on CL (continuous) - (eGFR/90)^theta
//   4) Race on VC - exp(theta_vc_{race}*I(race == {race}))

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  vector depot_1cmt_ode(real t, vector y, array[] real params, 
                        array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real ke = cl/vc;
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
   
    real slope = ka*y[1] - ke*y[2];
    real x = slope > 0 && y[2]/vc > y[4] ? slope/vc : 0;
    real z = t <= t_1 || (slope > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    vector[5] dydt;

    dydt[1] = -ka*y[1];                            // depot
    dydt[2] = slope;                               // central
    dydt[3] = t >= t_1 && t <= t_2 ? y[2] : 0;     // AUC_t_1-t_2
    dydt[4] = x;                                   // C_max
    dydt[5] = z;                                   // t_max
    
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
  
  array[n_subjects_new] real<lower = 0> wt;                    // baseline bodyweight (kg)
  array[n_subjects_new] int<lower = 0, upper = 1> cmppi;       // cmppi
  array[n_subjects_new] real<lower = 0> egfr;                  // eGFR
  int<lower = 2> n_races;                                      // number of unique races
  array[n_subjects_new] int<lower = 1, upper = n_races> race;  // race
  
  real<lower = 0> t_1;   // Time at which to start SS calculations (AUC_ss, C_max_ss, ...)
  real<lower = t_1> t_2; // Time at which to end SS calculations (AUC_ss, C_max_ss, ...)
  
  int<lower = 0, upper = 1> want_auc_cmax; // Want AUC and Cmax? If so, it'll 
                                           // use the ODE solution. Otherwise,
                                           // it'll use the analytical solution (and be faster)
  
}
transformed data{ 
  
  int n_random = 3; // Number of random effects
  int n_cmt = want_auc_cmax ? 5 : 2; // Number of compartments - depot, central (AUC_ss, Cmax_ss, Tmax_ss))
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt);
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  array[1, 2] real x_r = {{t_1, t_2}};

}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVKA;
  
  real theta_cl_wt;
  real theta_vc_wt;
  real theta_ka_cmppi;
  real theta_cl_egfr;
  real theta_vc_race2;
  real theta_vc_race3;
  real theta_vc_race4;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{

  vector[n_time_new] epred_stan; // f(TVs, x, eta = eta_new), eta_new ~ multi_normal(0, Omega) 
  vector[n_time_new] epred;      // epred_stan + error
  
  vector[want_auc_cmax ? n_subjects_new : 0] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[want_auc_cmax ? n_subjects_new : 0] c_max;   // Cmax between t1 and t2 (c_max_ss)
  vector[want_auc_cmax ? n_subjects_new : 0] t_max;   // Tmax between t1 and t2, then subtract off t1
  vector[n_subjects_new] t_half;                      // half-life
  
  vector[n_subjects_new] CL;
  vector[n_subjects_new] VC;
  vector[n_subjects_new] KA;
  vector[n_subjects_new] KE;

  {
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA});
    
    matrix[n_subjects_new, n_random] eta_new;
    matrix[n_subjects_new, n_random] theta_new;
    matrix[n_time_new, n_cmt] x_epred;

    for(i in 1:n_subjects_new){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects_new) .* exp(eta_new));
    
    vector[n_races] theta_vc_race = [0, theta_vc_race2, theta_vc_race3, 
                                     theta_vc_race4]';

    for(j in 1:n_subjects_new){
      
      real wt_over_70 = wt[j]/70;
      real wt_adjustment_cl = wt_over_70^theta_cl_wt;
      real wt_adjustment_vc = wt_over_70^theta_vc_wt;
      real cmppi_adjustment_ka = exp(theta_ka_cmppi*cmppi[j]);
      real egfr_adjustment_cl = (egfr[j]/90)^theta_cl_egfr;
      real race_adjustment_vc = exp(theta_vc_race[race[j]]);
      
      row_vector[n_random] theta_j_new = theta_new[j]; // access the parameters for subject j's epred
      
      CL[j] = theta_j_new[1];
      VC[j] = theta_j_new[2];
      KA[j] = theta_j_new[3];
      KE[j] = CL[j]/VC[j];
      
      CL[j] = theta_j_new[1] * wt_adjustment_cl * egfr_adjustment_cl;
      VC[j] = theta_j_new[2] * wt_adjustment_vc * race_adjustment_vc;
      KA[j] = theta_j_new[3] * cmppi_adjustment_ka;
      KE[j] = CL[j]/VC[j];
    
      if(want_auc_cmax == 1){
        
        x_epred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(depot_1cmt_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j], KA[j]}, bioav, tlag, x_r)';
                        
        epred_stan[subj_start[j]:subj_end[j]] = 
          x_epred[subj_start[j]:subj_end[j], 2] ./ VC[j];
          
        auc_ss[j] = max(x_epred[subj_start[j]:subj_end[j], 3]) / VC[j];
        c_max[j] = max(x_epred[subj_start[j]:subj_end[j], 4]);
        t_max[j] = max(x_epred[subj_start[j]:subj_end[j], 5]) - t_1;
          
      }else{
        
        x_epred[subj_start[j]:subj_end[j],] =
          pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           {CL[j], VC[j], KA[j]})';

        epred_stan[subj_start[j]:subj_end[j]] = 
          x_epred[subj_start[j]:subj_end[j], 2] ./ VC[j];
        
      }
  
      t_half[j] = log(2) ./ KE[j];
    
    }

    for(i in 1:n_time_new){
      if(epred_stan[i] == 0){
        epred[i] = 0;
      }else{
        real epred_tmp = epred_stan[i];
        real sigma_tmp_e = epred_tmp*sigma_p;
        epred[i] = normal_lb_rng(epred_tmp, sigma_tmp_e, 0.0);
      }
    }
  }
}


