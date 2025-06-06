// Some subjects have First Order Absorption (oral/subcutaneous), some IV, some 
//   both
// One-compartment PK Model
// IIV on CL, VC, KA, BIOAV (full covariance matrix)
// exponential error - DV = IPRED*exp(eps)
// User's choice of linear ODE solution or general ODE solution
// General ODE solution using Torsten will get out individual estimates of AUC, 
//   Cmax, Tmax, ... linear ODE will not

functions{
  
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
  
  int<lower = 0, upper = 1> want_auc_cmax; // Want AUC and Cmax? If so, it'll 
                                           // use the ODE solution. Otherwise,
                                           // it'll use the analytical solution (and be faster)
 
}
transformed data{ 
  
  int n_random = 4; // Number of random effects
  int n_cmt = want_auc_cmax ? 5 : 2; // Number of compartments - depot, central (AUC_ss, Cmax_ss, Tmax_ss))
  
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  array[1, 2] real x_r = {{t_1, t_2}};

}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = TVCL/TVVC> TVKA;
  real<lower = 0, upper = 1> TVBIOAV;
  
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
  
  vector[want_auc_cmax ? n_subjects : 0] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[want_auc_cmax ? n_subjects : 0] c_max;   // Cmax between t1 and t2 (c_max_ss)
  vector[want_auc_cmax ? n_subjects : 0] t_max;   // Tmax between t1 and t2, then subtract off t1
  vector[n_subjects] t_half;                      // half-life
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  vector[n_subjects] BIOAV;
  vector[n_subjects] KE;
 
  {
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA,
                                                         TVBIOAV});

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
    vector[n_subjects] KA_new;
    vector[n_subjects] BIOAV_new;
    vector[n_subjects] KE_new;

    matrix[n_time_new, 2] x_pred;
    matrix[n_time_new, 2] x_epred;
    matrix[n_time_new, want_auc_cmax ? 5 : 2] x_ipred;
    
    for(j in 1:n_subjects){
      
      row_vector[n_random] theta_j = theta[j]; // access the parameters for subject j
      row_vector[n_random] theta_j_new = theta_new[j]; // access the parameters for subject j's epred
      
      real cl_p = TVCL;
      real vc_p = TVVC;
      real ka_p = TVKA;
      real bioav_p = TVBIOAV; // Phi(inv_Phi(TVBIOAV) + covariate_effects);
      
      CL[j] = theta_j[1];
      VC[j] = theta_j[2];
      KA[j] = theta_j[3];
      BIOAV[j] = Phi(inv_Phi(TVBIOAV) + eta[j, 4]); // Phi(inv_Phi(TVBIOAV) + covariate_effects + eta[j, 4]);
      KE[j] = CL[j]/VC[j];
      
      CL_new[j] = theta_j_new[1];
      VC_new[j] = theta_j_new[2];
      KA_new[j] = theta_j_new[3];
      BIOAV_new[j] = Phi(inv_Phi(TVBIOAV) + eta_new[j, 4]); // Phi(inv_Phi(TVBIOAV) + covariate_effects + eta_new[j, 4]);
      KE_new[j] = CL_new[j]/VC_new[j];
      
      matrix[2, 2] K_epred = rep_matrix(0, 2, 2);
      matrix[2, 2] K_p = rep_matrix(0, 2, 2);
      
      K_p[1, 1] = -ka_p;
      K_p[2, 1] = ka_p;
      K_p[2, 2] = -cl_p/vc_p;

      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         K_p, {bioav_p, 1}, tlag[1:2])';
      
      K_epred[1, 1] = -KA_new[j];
      K_epred[2, 1] = KA_new[j];
      K_epred[2, 2] = -KE_new[j];
      
      x_epred[subj_start[j]:subj_end[j], ] =
        pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         K_epred, {BIOAV_new[j], 1}, tlag[1:2])';
      
      // x_pred[subj_start[j]:subj_end[j],] =
      //   pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
      //                    amt[subj_start[j]:subj_end[j]],
      //                    rate[subj_start[j]:subj_end[j]],
      //                    ii[subj_start[j]:subj_end[j]],
      //                    evid[subj_start[j]:subj_end[j]],
      //                    cmt[subj_start[j]:subj_end[j]],
      //                    addl[subj_start[j]:subj_end[j]],
      //                    ss[subj_start[j]:subj_end[j]],
      //                    {cl_p, vc_p, ka_p},
      //                    {bioav_p, 1})';
      // 
      // x_epred[subj_start[j]:subj_end[j],] =
      //   pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
      //                    amt[subj_start[j]:subj_end[j]],
      //                    rate[subj_start[j]:subj_end[j]],
      //                    ii[subj_start[j]:subj_end[j]],
      //                    evid[subj_start[j]:subj_end[j]],
      //                    cmt[subj_start[j]:subj_end[j]],
      //                    addl[subj_start[j]:subj_end[j]],
      //                    ss[subj_start[j]:subj_end[j]],
      //                    {CL_new[j], VC_new[j], KA_new[j]},
      //                    {BIOAV_new[j], 1})';
      
      pred[subj_start[j]:subj_end[j]] = 
        x_pred[subj_start[j]:subj_end[j], 2] ./ vc_p;

      epred_stan[subj_start[j]:subj_end[j]] = 
        x_epred[subj_start[j]:subj_end[j], 2] ./ VC_new[j];
      
      
      if(want_auc_cmax == 1){
        
        x_ipred[subj_start[j]:subj_end[j],] =
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
                         {CL[j], VC[j], KA[j]}, 
                         {BIOAV_new[j], 1, 1, 1, 1}, tlag, x_r)';
                        
        ipred[subj_start[j]:subj_end[j]] = 
          x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
          
        auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], 3]) / VC[j];
        c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 4]);
        t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]) - t_1;
          
      }else{
        
        matrix[2, 2] K = rep_matrix(0, 2, 2);
        
        K[1, 1] = -KA[j];
        K[2, 1] = KA[j];
        K[2, 2] = -CL[j]/VC[j];
      
        x_ipred[subj_start[j]:subj_end[j], ] =
          pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           K, {BIOAV[j], 1}, tlag[1:2])';
        
        // x_ipred[subj_start[j]:subj_end[j],] =
        //   pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
        //                    amt[subj_start[j]:subj_end[j]],
        //                    rate[subj_start[j]:subj_end[j]],
        //                    ii[subj_start[j]:subj_end[j]],
        //                    evid[subj_start[j]:subj_end[j]],
        //                    cmt[subj_start[j]:subj_end[j]],
        //                    addl[subj_start[j]:subj_end[j]],
        //                    ss[subj_start[j]:subj_end[j]],
        //                    {CL[j], VC[j], KA[j]},
        //                    {BIOAV[j], 1})';

        ipred[subj_start[j]:subj_end[j]] = 
          x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
        
      }

      t_half[j] = log(2) ./ KE[j];

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
