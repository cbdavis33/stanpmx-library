// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model with Indirect Response 4 (stimulation of dissipation) 
//   PD Model
// IIV on CL, VC, and KA (full covariance matrix)
// IIV on KIN, KOUT, SC50, and SMAX (full covariance matrix) for PD. HILL is 
//   fixed to be 1
// proportional plus additive error for PK - DV = IPRED*(1 + eps_p_pk) + eps_a_pk
// proportional plus additive error for PD - DV = IPRED*(1 + eps_p_pd) + eps_a_pd
// Observations are generated from a normal that is truncated below at 0
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   generates values from a normal that is truncated below at 0
// User's choice of coupled or non-coupled solver. Coupled solver doesn't output
//   Cmax, Tmax, ...
// PK output can include individual Cmax over the whole time period, Tmax between 
//   t1 and t2, AUC since 0 for every timepoint, and AUC between t1 and t2 (like 
//   a dosing interval)
// PD output can include individual Rmin (response decreases to a min), and Tmin
//   betweem t1 and t2 (like the time within a dosing interval at which the 
//   minimum is reached)

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  vector depot_1cmt_ir4_ode(real t, vector y, array[] real params, 
                            array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real kin = params[4];
    real kout = params[5];
    real sc50 = params[6];
    real smax = params[7];   
    real hill = params[8];   // It's fixed to 1 in this particular model
    real r_0 = params[9];
    
    real ke = cl/vc;
    
    real conc = y[2]/vc;
    
    real stim = (smax*pow(conc, hill))/(pow(sc50, hill) + pow(conc, hill));
    real response = y[3] + r_0;
    
    vector[3] dydt;

    dydt[1] = -ka*y[1];
    dydt[2] = ka*y[1] - ke*y[2];
    dydt[3] = kin - kout*(1 + stim)*response;
    
    return dydt;
    
  }
  
  vector depot_1cmt_ir4_ode_coupled(real t, vector y, vector y_pk, 
                                    array[] real params, array[] real x_r, 
                                    array[] int x_i){
    
    real vc = params[2];

    real kin = params[4];
    real kout = params[5];
    real sc50 = params[6];
    real smax = params[7];   
    real hill = params[8];   // It's fixed to 1 in this particular model
    real r_0 = params[9];
    
    real conc = y_pk[2]/vc;
    
    real stim = (smax*pow(conc, hill))/(pow(sc50, hill) + pow(conc, hill));
    real response = y[1] + r_0;
    
    vector[1] dydt;

    dydt[1] = kin - kout*(1 + stim)*response;
    
    return dydt;
    
  }
  

  vector depot_1cmt_ir4_with_auc_cmax_ode(real t, vector y, array[] real params, 
                                          array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real kin = params[4];
    real kout = params[5];
    real sc50 = params[6];
    real smax = params[7];   
    real hill = params[8];   // It's fixed to 1 in this particular model
    real r_0 = params[9];
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real ke = cl/vc;
    
    real conc = y[2]/vc;
    
    real stim = (smax*pow(conc, hill))/(pow(sc50, hill) + pow(conc, hill));
    real response = y[3] + r_0;
    
    real slope_pk = ka*y[1] - ke*y[2];
    real x_pk = slope_pk > 0 && conc > y[6] ? slope_pk/vc : 0;
    real z_pk = t <= t_1 || (slope_pk > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    real slope_pd = kin - kout*(1 + stim)*response;
    real x_pd = slope_pd < 0 && y[3] < y[8] ? slope_pd : 0;
    real z_pd = t <= t_1 || (slope_pd < 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    vector[9] dydt;

    dydt[1] = -ka*y[1];
    dydt[2] = slope_pk;
    dydt[3] = slope_pd;
    dydt[4] = y[2];                                // AUC
    dydt[5] = t >= t_1 && t <= t_2 ? y[2] : 0;     // AUC_t_1-t_2
    dydt[6] = x_pk;                                // C_max 
    dydt[7] = z_pk;                                // t_max for PK
    dydt[8] = x_pd;                                // R_max
    dydt[9] = z_pd;                                // t_max for PD
    
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
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_ka;
  
  real<lower = 0> TVKIN;
  real<lower = 0> TVKOUT;
  real<lower = 0> TVSC50;
  real<lower = 0> TVSMAX;
  
  real<lower = 0> omega_kin;
  real<lower = 0> omega_kout;
  real<lower = 0> omega_sc50;
  real<lower = 0> omega_smax;
  
  corr_matrix[3] R_pk;  // Correlation matrix before transforming to Omega.
                        // Can in theory change this to having inputs for
                        // cor_cl_vc, cor_cl_ka, ... and then construct the 
                        // correlation matrix in transformed data, but it's easy
                        // enough to do in R
                        
  corr_matrix[4] R_pd;  // Correlation matrix before transforming to Omega.
                        // Can in theory change this to having inputs for
                        // cor_kin_kout, cor_kin_ic50, ... and then construct the 
                        // correlation matrix in transformed data, but it's easy
                        // enough to do in R
  
  real<lower = 0> sigma_p_pk;
  real<lower = 0> sigma_a_pk;
  real<lower = -1, upper = 1> cor_p_a_pk;
  
  real<lower = 0> sigma_p_pd;
  real<lower = 0> sigma_a_pd;
  real<lower = -1, upper = 1> cor_p_a_pd;
  
  int<lower = 0, upper = 1> coupled;
  int<lower = 0, upper = 1> want_auc_cmax; // non-coupled only
  int<lower = 0, upper = 3> solver;        // 1 = rk45, 2 = bdf, 3 = adams (non-coupled only)
  
  real<lower = 0> t_1;
  real<lower = t_1> t_2;
   
}
transformed data{
  
  int n_random_pk = 3;
  int n_random_pd = 4;
  
  int n_cmt_pk = 2;    // number of ODEs in PK model (depot, central)
  int n_cmt_pd = 1;    // number of ODEs in PD system
  int n_cmt_extra = 6; // number of ODEs for AUC, Tmax, Rmin, ...
  int n_cmt = want_auc_cmax ? n_cmt_pk + n_cmt_pd + n_cmt_extra : n_cmt_pk + n_cmt_pd;

  vector[n_random_pk] omega_pk = [omega_cl, omega_vc, omega_ka]';
                               
  vector[n_random_pd] omega_pd = [omega_kin, omega_kout, omega_sc50, omega_smax]';
  
  matrix[n_random_pk, n_random_pk] L_pk = cholesky_decompose(R_pk);
  
  vector[2] sigma_pk = [sigma_p_pk, sigma_a_pk]';
  matrix[2, 2] R_Sigma_pk = rep_matrix(1, 2, 2);
  R_Sigma_pk[1, 2] = cor_p_a_pk;
  R_Sigma_pk[2, 1] = cor_p_a_pk;
  
  matrix[2, 2] Sigma_pk = quad_form_diag(R_Sigma_pk, sigma_pk);
  
  matrix[n_random_pd, n_random_pd] L_pd = cholesky_decompose(R_pd);
  
  vector[2] sigma_pd = [sigma_p_pd, sigma_a_pd]';
  matrix[2, 2] R_Sigma_pd = rep_matrix(1, 2, 2);
  R_Sigma_pd[1, 2] = cor_p_a_pd;
  R_Sigma_pd[2, 1] = cor_p_a_pd;
  
  matrix[2, 2] Sigma_pd = quad_form_diag(R_Sigma_pd, sigma_pd);
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);  // Hardcoding, but could be data or a parameter in another situation
  
  real TVHILL = 1.0;                                                   // HILL is fixed to 1.0 in this model,
  vector<lower = 0>[n_subjects] HILL = rep_vector(TVHILL, n_subjects); // but it could be data or a parameter in another model.
                                                                       // Putting this here will require the least amount of 
                                                                       // changes in the code if that change is made
  
  array[1, 2] real x_r = {{t_1, t_2}};
  
}
model{
  
}
generated quantities{
  
  vector[n_total] ipred; // concentration or response with no residual error
  vector[n_total] dv;    // concentration or response with residual error
  
  vector[want_auc_cmax ? n_total : 0] auc;            // AUC 
  vector[want_auc_cmax ? n_subjects : 0] auc_t1_t2;   // AUC from t1 up to t2
  vector[want_auc_cmax ? n_subjects : 0] c_max;       // Cmax
  vector[want_auc_cmax ? n_subjects : 0] t_max;       // Tmax foor PK
  vector[want_auc_cmax ? n_subjects : 0] t_half;      // half-life
  vector[want_auc_cmax ? n_subjects : 0] r_min;       // Minimum response
  vector[want_auc_cmax ? n_subjects : 0] t_min_pd;    // Tmin for PD
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  
  vector[n_subjects] KIN;
  vector[n_subjects] KOUT;
  vector[n_subjects] SC50;
  vector[n_subjects] SMAX;
  
  {
  
    row_vector[n_random_pk] typical_values_pk = to_row_vector({TVCL, TVVC, TVKA});
    
    matrix[n_subjects, n_random_pk] eta_pk;   
    matrix[n_subjects, n_random_pk] theta_pk; 
    
    row_vector[n_random_pd] typical_values_pd = to_row_vector({TVKIN, TVKOUT, 
                                                               TVSC50, TVSMAX});
    
    matrix[n_subjects, n_random_pd] eta_pd;   
    matrix[n_subjects, n_random_pd] theta_pd; 
  
    matrix[n_total, n_cmt] x_ipred;
    
    vector[n_subjects] r_0;
    
    for(i in 1:n_subjects){
      eta_pk[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random_pk),
                                              diag_pre_multiply(omega_pk, L_pk))';
                                           
      eta_pd[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random_pd),
                                              diag_pre_multiply(omega_pd, L_pd))';                                           
    }
    theta_pk = (rep_matrix(typical_values_pk, n_subjects) .* exp(eta_pk));
    theta_pd = (rep_matrix(typical_values_pd, n_subjects) .* exp(eta_pd));

    CL = col(theta_pk, 1);
    VC = col(theta_pk, 2);
    KA = col(theta_pk, 3);
    
    KIN = col(theta_pd, 1);
    KOUT = col(theta_pd, 2);
    SC50 = col(theta_pd, 3);
    SMAX = col(theta_pd, 4);
    
    r_0 = KIN ./ KOUT; 
    
    for(j in 1:n_subjects){
      if(coupled){
        if(want_auc_cmax){
          reject("This doesn't return AUC, Cmax, ... with the coupled solver.");
        }else{
          if(solver == 1){
            
            x_ipred[subj_start[j]:subj_end[j],] =
              pmx_solve_onecpt_rk45(depot_1cmt_ir4_ode_coupled,
                                    n_cmt_pd,
                                    time[subj_start[j]:subj_end[j]],
                                    amt[subj_start[j]:subj_end[j]],
                                    rate[subj_start[j]:subj_end[j]],
                                    ii[subj_start[j]:subj_end[j]],
                                    evid[subj_start[j]:subj_end[j]],
                                    cmt[subj_start[j]:subj_end[j]],
                                    addl[subj_start[j]:subj_end[j]],
                                    ss[subj_start[j]:subj_end[j]],
                                    {CL[j], VC[j], KA[j], 
                                     KIN[j], KOUT[j], SC50[j], SMAX[j], 
                                     HILL[j], r_0[j]}, 
                                    bioav, tlag)';
                                    
          }else if(solver == 2){
            
            x_ipred[subj_start[j]:subj_end[j],] =
              pmx_solve_onecpt_bdf(depot_1cmt_ir4_ode_coupled,
                                   n_cmt_pd,
                                   time[subj_start[j]:subj_end[j]],
                                   amt[subj_start[j]:subj_end[j]],
                                   rate[subj_start[j]:subj_end[j]],
                                   ii[subj_start[j]:subj_end[j]],
                                   evid[subj_start[j]:subj_end[j]],
                                   cmt[subj_start[j]:subj_end[j]],
                                   addl[subj_start[j]:subj_end[j]],
                                   ss[subj_start[j]:subj_end[j]],
                                   {CL[j], VC[j], KA[j], 
                                    KIN[j], KOUT[j], SC50[j], SMAX[j], 
                                    HILL[j], r_0[j]}, 
                                   bioav, tlag)';
            
          }else{
            reject("The adams solver isn't available for the coupled solver.");
          }
        }
      }else{
        if(want_auc_cmax){
          if(solver == 1){
            
            x_ipred[subj_start[j]:subj_end[j],] =
              pmx_solve_rk45(depot_1cmt_ir4_with_auc_cmax_ode,
                             n_cmt,
                             time[subj_start[j]:subj_end[j]],
                             amt[subj_start[j]:subj_end[j]],
                             rate[subj_start[j]:subj_end[j]],
                             ii[subj_start[j]:subj_end[j]],
                             evid[subj_start[j]:subj_end[j]],
                             cmt[subj_start[j]:subj_end[j]],
                             addl[subj_start[j]:subj_end[j]],
                             ss[subj_start[j]:subj_end[j]],
                             {CL[j], VC[j], KA[j], 
                              KIN[j], KOUT[j], SC50[j], SMAX[j], 
                              HILL[j], r_0[j]},
                             bioav, tlag, x_r)';
                             
          }else if(solver == 2){
            
            x_ipred[subj_start[j]:subj_end[j],] =
              pmx_solve_bdf(depot_1cmt_ir4_with_auc_cmax_ode,
                             n_cmt,
                             time[subj_start[j]:subj_end[j]],
                             amt[subj_start[j]:subj_end[j]],
                             rate[subj_start[j]:subj_end[j]],
                             ii[subj_start[j]:subj_end[j]],
                             evid[subj_start[j]:subj_end[j]],
                             cmt[subj_start[j]:subj_end[j]],
                             addl[subj_start[j]:subj_end[j]],
                             ss[subj_start[j]:subj_end[j]],
                             {CL[j], VC[j], KA[j], 
                              KIN[j], KOUT[j], SC50[j], SMAX[j], 
                              HILL[j], r_0[j]},
                             bioav, tlag, x_r)';
            
          }else{
            
            x_ipred[subj_start[j]:subj_end[j],] =
              pmx_solve_adams(depot_1cmt_ir4_with_auc_cmax_ode,
                              n_cmt,
                              time[subj_start[j]:subj_end[j]],
                              amt[subj_start[j]:subj_end[j]],
                              rate[subj_start[j]:subj_end[j]],
                              ii[subj_start[j]:subj_end[j]],
                              evid[subj_start[j]:subj_end[j]],
                              cmt[subj_start[j]:subj_end[j]],
                              addl[subj_start[j]:subj_end[j]],
                              ss[subj_start[j]:subj_end[j]],
                              {CL[j], VC[j], KA[j], 
                               KIN[j], KOUT[j], SC50[j], SMAX[j], 
                               HILL[j], r_0[j]},
                              bioav, tlag, x_r)';
            
          }
          
          auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 4] ./ VC[j];
          
          auc_t1_t2[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]) / VC[j];
          c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 6]);
          t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 7]) - t_1;
          t_half[j] = log(2)/(CL[j]/VC[j]);
      
          r_min[j] = min(x_ipred[subj_start[j]:subj_end[j], 8]) + r_0[j];
          t_min_pd[j] = max(x_ipred[subj_start[j]:subj_end[j], 9]) - t_1;
          
        }else{
          if(solver == 1){
            
            x_ipred[subj_start[j]:subj_end[j],] =
              pmx_solve_rk45(depot_1cmt_ir4_ode,
                             n_cmt,
                             time[subj_start[j]:subj_end[j]],
                             amt[subj_start[j]:subj_end[j]],
                             rate[subj_start[j]:subj_end[j]],
                             ii[subj_start[j]:subj_end[j]],
                             evid[subj_start[j]:subj_end[j]],
                             cmt[subj_start[j]:subj_end[j]],
                             addl[subj_start[j]:subj_end[j]],
                             ss[subj_start[j]:subj_end[j]],
                             {CL[j], VC[j], KA[j], 
                              KIN[j], KOUT[j], SC50[j], SMAX[j], 
                              HILL[j], r_0[j]},
                             bioav, tlag, x_r)';
            
          }else if(solver == 2){
            
            x_ipred[subj_start[j]:subj_end[j],] =
              pmx_solve_bdf(depot_1cmt_ir4_ode,
                            n_cmt,
                            time[subj_start[j]:subj_end[j]],
                            amt[subj_start[j]:subj_end[j]],
                            rate[subj_start[j]:subj_end[j]],
                            ii[subj_start[j]:subj_end[j]],
                            evid[subj_start[j]:subj_end[j]],
                            cmt[subj_start[j]:subj_end[j]],
                            addl[subj_start[j]:subj_end[j]],
                            ss[subj_start[j]:subj_end[j]],
                            {CL[j], VC[j], KA[j], 
                             KIN[j], KOUT[j], SC50[j], SMAX[j], 
                             HILL[j], r_0[j]},
                            bioav, tlag, x_r)';
            
          }else{
            
            x_ipred[subj_start[j]:subj_end[j],] =
              pmx_solve_adams(depot_1cmt_ir4_ode,
                              n_cmt,
                              time[subj_start[j]:subj_end[j]],
                              amt[subj_start[j]:subj_end[j]],
                              rate[subj_start[j]:subj_end[j]],
                              ii[subj_start[j]:subj_end[j]],
                              evid[subj_start[j]:subj_end[j]],
                              cmt[subj_start[j]:subj_end[j]],
                              addl[subj_start[j]:subj_end[j]],
                              ss[subj_start[j]:subj_end[j]],
                              {CL[j], VC[j], KA[j], 
                               KIN[j], KOUT[j], SC50[j], SMAX[j], 
                               HILL[j], r_0[j]},
                              bioav, tlag, x_r)';
            
          }
        }
      }

      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){
          ipred[k] = x_ipred[k, 2] / VC[j];
        }else if(cmt[k] == 3){
          ipred[k] = x_ipred[k, 3] + r_0[j];
        }
      }

    }
    

    for(i in 1:n_total){
      if(ipred[i] == 0){
         dv[i] = 0;
      }else{
        
        if(cmt[i] == 2 || cmt[i] == 3){
          real ipred_tmp = ipred[i];
          real sigma_tmp = cmt[i] == 2 ? 
            sqrt(square(ipred_tmp) * Sigma_pk[1, 1] + Sigma_pk[2, 2] + 
                 2*ipred_tmp*Sigma_pk[2, 1]) :
            sqrt(square(ipred_tmp) * Sigma_pd[1, 1] + Sigma_pd[2, 2] + 
                 2*ipred_tmp*Sigma_pd[2, 1]);
          
          dv[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
          
        }
      }
    }
  }
}

