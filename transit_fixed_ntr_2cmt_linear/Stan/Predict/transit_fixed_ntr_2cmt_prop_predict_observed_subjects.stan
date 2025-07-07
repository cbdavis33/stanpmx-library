// Two-compartment PK Model with Transit Compartment absorption
// IIV on CL, VC, Q, VP, KA, MTT (full covariance matrix)
// n_transit is fixed to a positive integer that can be defined as data
// The 0th transit compartment is where the dosing happens (cmt = 1). This could
//   also be called the Depot
// User's choice of linear ODE solution or general ODE solution
// General ODE solution using Torsten will get out individual estimates of AUC, 
//   Cmax, Tmax, ... linear ODE will not
// proportional error - DV = IPRED*(1 + eps_p)
// Predictions are generated from a normal that is truncated below at 0

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  vector transit_fixed_2cmt_ode(real t, vector y, array[] real params, 
                                array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real q = params[3];
    real vp = params[4];
    real ka = params[5];
    real ktr = params[6];
    
    int n_transit = x_i[1];
    
    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;
    
    vector[n_transit + 4] dydt;

    dydt[1] = -ktr*y[1];               // 0th transit compartment (depot)
    for(i in 2:(n_transit + 1)){
      dydt[i] = ktr*(y[i - 1] - y[i]);
    }
    dydt[n_transit + 2] = ktr*y[n_transit + 1] - ka*y[n_transit + 2];          // absorption
    dydt[n_transit + 3] = ka*y[n_transit + 2] - (ke + k_cp)*y[n_transit + 3] + // central
                          k_pc*y[n_transit + 4];
    dydt[n_transit + 4] = k_cp*y[n_transit + 3] - k_pc*y[n_transit + 4];       // peripheral
    
    return dydt;
  }
  
  vector transit_fixed_2cmt_with_auc_ode(real t, vector y, array[] real params,
                                         array[] real x_r, array[] int x_i){

    real cl = params[1];
    real vc = params[2];
    real q = params[3];
    real vp = params[4];
    real ka = params[5];
    real ktr = params[6];
    
    int n_transit = x_i[1];
    
    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;

    real t_1 = x_r[1];
    real t_2 = x_r[2];

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
  
  int<lower = 0, upper = 1> want_auc_cmax; // Want AUC and Cmax? If so, it'll 
                                           // use the ODE solution. Otherwise,
                                           // it'll use the linear ODE solution (and be faster)
 
}
transformed data{ 
  
  int n_random = 6;           // Number of random effects
  int n_cmt = want_auc_cmax ? n_transit + 8 : n_transit + 4; // Depot, tr_1, ..., tr_n, absorption, central, peripheral
                                                             //   (AUC, AUC_t1-t2, Cmax, Tmax)

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
  real<lower = TVCL/TVVC> TVKA;
  real<lower = 0> TVMTT;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{
  
  vector[n_time_new] pred;       // f(TVs, x, eta = 0) 
  vector[n_time_new] epred_stan; // f(TVs, x, eta = eta_new), eta_new ~ multi_normal(0, Omega) 
  vector[n_time_new] epred;      // epred_stan + error
  vector[n_time_new] ipred;      // f(TVs, x, eta = eta_i), eta_i are etas for observed subjects
  vector[n_time_new] dv;         // ipred + error
  
  vector[want_auc_cmax ? n_time_new : 0] auc;     // AUC
  vector[want_auc_cmax ? n_subjects : 0] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[want_auc_cmax ? n_subjects : 0] c_max;   // Cmax between t1 and t2 (c_max_ss)
  vector[want_auc_cmax ? n_subjects : 0] t_max;   // Tmax between t1 and t2, then subtract off t1
  vector[n_subjects] t_half_alpha;                // alpha half-life
  vector[n_subjects] t_half_terminal;             // terminal half-life
  
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
    vector[n_subjects] Q_new;
    vector[n_subjects] VP_new;
    vector[n_subjects] KA_new;
    vector[n_subjects] MTT_new;
    vector[n_subjects] KTR_new;

    matrix[n_time_new, n_transit + 4] x_pred;
    matrix[n_time_new, n_transit + 4] x_epred;
    matrix[n_time_new, want_auc_cmax ? n_transit + 8 : n_transit + 4] x_ipred;
    
    for(j in 1:n_subjects){
      
      row_vector[n_random] theta_j = theta[j]; // access the parameters for subject j
      row_vector[n_random] theta_j_new = theta_new[j]; // access the parameters for subject j's epred
      
      real cl_p = TVCL;
      real vc_p = TVVC;
      real q_p = TVQ;
      real vp_p = TVVP;
      real ka_p = TVKA;
      real mtt_p = TVMTT;
      real ktr_p = (n_transit + 1) ./ mtt_p;
      
      CL[j] = theta_j[1];
      VC[j] = theta_j[2];
      Q[j] = theta_j[3];
      VP[j] = theta_j[4];
      KA[j] = theta_j[5];
      MTT[j] = theta_j[6];
      KTR[j] = (n_transit + 1)/MTT[j];

      CL_new[j] = theta_j_new[1];
      VC_new[j] = theta_j_new[2];
      Q_new[j] = theta_j_new[3];
      VP_new[j] = theta_j_new[4];
      KA_new[j] = theta_j_new[5];
      MTT_new[j] = theta_j_new[6];
      KTR_new[j] = (n_transit + 1)/MTT_new[j];
      
      matrix[n_transit + 4, n_transit + 4] K_epred = rep_matrix(0, n_transit + 4, 
                                                                n_transit + 4);
      matrix[n_transit + 4, n_transit + 4] K_p = rep_matrix(0, n_transit + 4, 
                                                            n_transit + 4);
                      
      for(i in 1:(n_transit + 1)){
        K_epred[i, i] = -KTR_new[j];
        K_epred[(i + 1), i] = KTR_new[j];
      }
        
      K_epred[(n_transit + 2), (n_transit + 2)] = -KA_new[j];
      K_epred[(n_transit + 3), (n_transit + 2)] = KA_new[j];
      K_epred[(n_transit + 3), (n_transit + 3)] = -((CL_new[j] + Q_new[j])/VC_new[j]);
      K_epred[(n_transit + 3), (n_transit + 4)] = Q_new[j]/VP_new[j];
      K_epred[(n_transit + 4), (n_transit + 3)] = Q_new[j]/VC_new[j];
      K_epred[(n_transit + 4), (n_transit + 4)] = -Q_new[j]/VP_new[j];
      
      x_epred[subj_start[j]:subj_end[j],] =
        pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         K_epred, bioav, tlag)';
        
      for(i in 1:(n_transit + 1)){
        K_p[i, i] = -ktr_p;
        K_p[(i + 1), i] = ktr_p;
      }
        
      K_p[(n_transit + 2), (n_transit + 2)] = -ka_p;
      K_p[(n_transit + 3), (n_transit + 2)] = ka_p;
      K_p[(n_transit + 3), (n_transit + 3)] = -((cl_p + q_p)/vc_p);
      K_p[(n_transit + 3), (n_transit + 4)] = q_p/vp_p;
      K_p[(n_transit + 4), (n_transit + 3)] = q_p/vc_p;
      K_p[(n_transit + 4), (n_transit + 4)] = -q_p/vp_p;
      
      x_pred[subj_start[j]:subj_end[j],] =
        pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         K_p, bioav, tlag)';
                         
      pred[subj_start[j]:subj_end[j]] = 
        x_pred[subj_start[j]:subj_end[j], n_transit + 3] ./ vc_p;

      epred_stan[subj_start[j]:subj_end[j]] = 
        x_epred[subj_start[j]:subj_end[j], n_transit + 3] ./ VC_new[j];
      
      if(want_auc_cmax == 1){
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(transit_fixed_2cmt_with_auc_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j], Q[j], VP[j], KA[j], KTR[j]},
                         bioav, tlag, x_r, x_i)';
        
        auc[subj_start[j]:subj_end[j]] = 
                    x_ipred[subj_start[j]:subj_end[j], (n_transit + 5)] ./ VC[j];
                    
        auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], n_transit + 6]) / VC[j];
        c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], n_transit + 7]);
        t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], n_transit + 8]) - t_1;
        
      }else{
        
        matrix[n_transit + 4, n_transit + 4] K = rep_matrix(0, n_transit + 4, 
                                                            n_transit + 4);
                                                            
        for(i in 1:(n_transit + 1)){
          K[i, i] = -KTR[j];
          K[(i + 1), i] = KTR[j];
        }
          
        K[(n_transit + 2), (n_transit + 2)] = -KA[j];
        K[(n_transit + 3), (n_transit + 2)] = KA[j];
        K[(n_transit + 3), (n_transit + 3)] = -((CL[j] + Q[j])/VC[j]);
        K[(n_transit + 3), (n_transit + 4)] = Q[j]/VP[j];
        K[(n_transit + 4), (n_transit + 3)] = Q[j]/VC[j];
        K[(n_transit + 4), (n_transit + 4)] = -Q[j]/VP[j];
      
        x_ipred[subj_start[j]:subj_end[j],] =
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
      
      ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], n_transit + 3] ./ VC[j];
      
    }
    
    vector[n_subjects] alpha = 0.5*(CL./VC + Q./VC + Q./VP + 
                            sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
    vector[n_subjects] beta = 0.5*(CL./VC + Q./VC + Q./VP - 
                            sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
    
    t_half_alpha = log(2)/alpha;
    t_half_terminal = log(2)/beta;

    for(i in 1:n_time_new){
      
      if(epred_stan[i] == 0){
        epred[i] = 0;
      }else{
        real epred_tmp = epred_stan[i];
        real sigma_tmp_e = epred_tmp*sigma_p;
        epred[i] = normal_lb_rng(epred_tmp, sigma_tmp_e, 0.0);
      }
      
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

