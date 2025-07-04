// One-compartment PK Model with Transit Compartment absorption
// IIV on CL, VC, KA, MTT (full covariance matrix)
// n_transit is fixed to a positive integer that can be defined as data
// The 0th transit compartment is where the dosing happens (cmt = 1). This could
//   also be called the Depot
// User's choice of linear ODE solution or general ODE solution
// General ODE solution using Torsten will get out individual estimates of AUC, 
//   Cmax, Tmax, ... linear ODE will not
// proportional plus additive error - DV = IPRED*(1 + eps_p) + eps_a
// Predictions are generated from a normal that is truncated below at 0

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  vector transit_fixed_1cmt_ode(real t, vector y, array[] real params, 
                                array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real ktr = params[4];
    
    int n_transit = x_i[1];
    
    real ke = cl/vc;
    
    vector[n_transit + 3] dydt;

    dydt[1] = -ktr*y[1];               // 0th transit compartment (depot)
    for(i in 2:(n_transit + 1)){
      dydt[i] = ktr*(y[i - 1] - y[i]);
    }
    dydt[n_transit + 2] = ktr*y[n_transit + 1] - ka*y[n_transit + 2]; // absorption
    dydt[n_transit + 3] = ka*y[n_transit + 2] - ke*y[n_transit + 3];  // central
    
    return dydt;
    
  }
  
  vector transit_fixed_1cmt_with_auc_ode(real t, vector y, array[] real params,
                                         array[] real x_r, array[] int x_i){

    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real ktr = params[4];

    real t_1 = x_r[1];
    real t_2 = x_r[2];

    int n_transit = x_i[1];

    real ke = cl/vc;

    real slope = ka*y[n_transit + 2] - ke*y[n_transit + 3];
    real x = slope > 0 && y[n_transit + 3]/vc > y[n_transit + 6] ? slope/vc : 0;
    real z = t <= t_1 || (slope > 0 && t >= t_1 && t <= t_2) ? 1 : 0;

    vector[n_transit + 7] dydt;

    dydt[1] = -ktr*y[1];                               // 0th transit compartment (depot)
    for(i in 2:(n_transit + 1)){
      dydt[i] = ktr*(y[i - 1] - y[i]);
    }
    dydt[n_transit + 2] = ktr*y[n_transit + 1] - ka*y[n_transit + 2];
    dydt[n_transit + 3] = slope;
    dydt[n_transit + 4] = y[n_transit + 3];                              // AUC
    dydt[n_transit + 5] = t >= t_1 && t <= t_2 ? y[n_transit + 3] : 0;   // AUC_t_1-t_2
    dydt[n_transit + 6] = x;                                             // C_max
    dydt[n_transit + 7] = z;                                             // t_max

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
  
  int<lower = 1> n_transit;
  
  real<lower = 0> t_1;   // Time at which to start SS calculations (AUC_ss, C_max_ss, ...)
  real<lower = t_1> t_2; // Time at which to end SS calculations (AUC_ss, C_max_ss, ...)
 
  int<lower = 0, upper = 1> want_auc_cmax; // Want AUC and Cmax? If so, it'll 
                                           // use the ODE solution. Otherwise,
                                           // it'll use the linear ODE solution (and be faster)
 
}
transformed data{ 
  
  int n_random = 4;           // Number of random effects
  int n_cmt = want_auc_cmax ? n_transit + 7 : n_transit + 3; // Depot, tr_1, ..., tr_n, absorption, central,
                                                             //   (AUC, AUC_t1-t2, Cmax, Tmax)
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt);
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  array[1, 2] real x_r = {{t_1, t_2}};
  
  array[1, 1] int x_i = {{n_transit}};
  
}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = TVCL/TVVC> TVKA;
  real<lower = 0> TVMTT;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  vector<lower = 0>[2] sigma;
  cholesky_factor_corr[2] L_Sigma;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{
  
  vector[n_time_new] epred_stan; // f(TVs, x, eta = eta_new), eta_new ~ multi_normal(0, Omega) 
  vector[n_time_new] epred;      // epred_stan + error
  
  vector[want_auc_cmax ? n_time_new : 0] auc;         // AUC from 0 up until t
  vector[want_auc_cmax ? n_subjects_new : 0] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[want_auc_cmax ? n_subjects_new : 0] c_max;   // Cmax between t1 and t2 (c_max_ss)
  vector[want_auc_cmax ? n_subjects_new : 0] t_max;   // Tmax between t1 and t2, then subtract off t1
  vector[n_subjects_new] t_half;                      // half-life
  
  vector[n_subjects_new] CL;
  vector[n_subjects_new] VC;
  vector[n_subjects_new] KA;
  vector[n_subjects_new] MTT;
  vector[n_subjects_new] KTR;
 
  {
    
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA, 
                                                         TVMTT});
                                                         
    matrix[n_subjects_new, n_random] eta_new;
    matrix[n_subjects_new, n_random] theta_new;
    matrix[n_time_new, n_cmt] x_epred;
    
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
    
    real sigma_sq_p = Sigma[1, 1];
    real sigma_sq_a = Sigma[2, 2];
    real sigma_p_a = Sigma[1, 2];

    for(i in 1:n_subjects_new){
      eta_new[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                               diag_pre_multiply(omega, L))';
    }
    theta_new = (rep_matrix(typical_values, n_subjects_new) .* exp(eta_new));
    
    for(j in 1:n_subjects_new){
      
      row_vector[n_random] theta_j_new = theta_new[j]; // access the parameters for subject j's epred
      
      CL[j] = theta_j_new[1];
      VC[j] = theta_j_new[2];
      KA[j] = theta_j_new[3];
      MTT[j] = theta_j_new[4];
      KTR[j] = (n_transit + 1)/MTT[j];
       
      t_half[j] = log(2)/(CL[j]/VC[j]);
      
      if(want_auc_cmax == 1){
        
        x_epred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(transit_fixed_1cmt_with_auc_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j], KA[j], KTR[j]}, bioav, tlag, x_r, x_i)';
          
        auc[subj_start[j]:subj_end[j]] = 
                    x_epred[subj_start[j]:subj_end[j], (n_transit + 4)] ./ VC[j];
                    
        auc_ss[j] = max(x_epred[subj_start[j]:subj_end[j], n_transit + 5]) / VC[j];
        c_max[j] = max(x_epred[subj_start[j]:subj_end[j], n_transit + 6]);
        t_max[j] = max(x_epred[subj_start[j]:subj_end[j], n_transit + 7]) - t_1;
        
      }else{
       
        matrix[n_transit + 3, n_transit + 3] K = rep_matrix(0, n_transit + 3, 
                                                            n_transit + 3);
                                                            
        for(i in 1:(n_transit + 1)){
          K[i, i] = -KTR[j];
          K[(i + 1), i] = KTR[j];
        }
          
        K[(n_transit + 2), (n_transit + 2)] = -KA[j];
        K[(n_transit + 3), (n_transit + 2)] = KA[j];
        K[(n_transit + 3), (n_transit + 3)] = -CL[j]/VC[j];
      
        x_epred[subj_start[j]:subj_end[j], ] =
          pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           K, bioav, tlag)';
       
      }
     
      epred_stan[subj_start[j]:subj_end[j]] = 
          x_epred[subj_start[j]:subj_end[j], n_transit + 3] ./ VC[j];
          
    }

    for(i in 1:n_time_new){
      if(epred_stan[i] == 0){
        epred[i] = 0;
      }else{
        real epred_tmp = epred_stan[i];
        real sigma_tmp_e = sqrt(square(epred_tmp) * sigma_sq_p + sigma_sq_a + 
                                2*epred_tmp*sigma_p_a);
        epred[i] = normal_lb_rng(epred_tmp, sigma_tmp_e, 0.0);
      }
    }
  }
}
