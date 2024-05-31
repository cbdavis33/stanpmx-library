// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model with Friberg-Karlsson Neutropenia Model
// IIV on CL, VC, and KA (full covariance matrix)
// IIV on MTT, CIRC0, GAMMA, ALPHA (full covariance matrix) for PD
// proportional plus additive error for PK - DV = IPRED*(1 + eps_p) + eps_a
// exponential error for PD - DV = IPRED*exp(eps_pd)
// General ODE solution using Torsten
// Since we have a normal distribution on the PK error, but the DV must be > 0, 
//   it generates values from a normal that is truncated below at 0
// PK output includes individual Cmax over the whole time period, Tmax between 
//   t1 and t2, AUC since 0 for every timepoint, and AUC between t1 and t2 (like 
//   a dosing interval)
// PD output includes individual Rmin (neutrophils decrease to a min)

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  vector depot_1cmt_friberg_ode(real t, vector y, array[] real params, 
                                array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real mtt = params[4];
    real circ_0 = params[5];
    real gamma = params[6];
    real alpha = params[7];
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real ke = cl/vc;
    real k_tr = 4/mtt; // k_tr = (n_tr + 1)/mtt    
    
    real conc = y[2]/vc;
    
    real e_drug = fmin(1.0, alpha*conc); // Maybe reparameterize this so no more fmin?
    real prol = y[3] + circ_0;
    real transit_1 = y[4] + circ_0; 
    real transit_2 = y[5] + circ_0;
    real transit_3 = y[6] + circ_0;
    real circ = y[7] + circ_0; // fmax(machine_precision(), y[7] + circ_0)
    
    real slope_pk = ka*y[1] - ke*y[2];
    real x_pk = slope_pk > 0 && conc > y[10] ? slope_pk/vc : 0;
    real z_pk = t <= t_1 || (slope_pk > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    real slope_pd = k_tr*(transit_3 - circ);
    real x_pd = slope_pd < 0 && y[7] < y[12] ? slope_pd : 0;
    
    vector[12] dydt;
    
    dydt[1] = -ka*y[1];                               // depot
    dydt[2] = slope_pk;                               // central
    dydt[3] = k_tr*prol*((1 - e_drug)*(circ_0/circ)^gamma - 1);  // proliferative cells
    dydt[4] = k_tr*(prol - transit_1);                // transit 1
    dydt[5] = k_tr*(transit_1 - transit_2);           // transit 2
    dydt[6] = k_tr*(transit_2 - transit_3);           // transit 3
    dydt[7] = slope_pd;                               // circulating blood cells
    dydt[8] = y[2];                                   // AUC
    dydt[9] = t >= t_1 && t <= t_2 ? y[2] : 0;        // AUC_t_1-t_2
    dydt[10] = x_pk;                                  // C_max 
    dydt[11] = z_pk;                                  // t_max for PK
    dydt[12] = x_pd;                                  // R_min
    
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
  
  real<lower = 0> TVMTT;
  real<lower = 0> TVCIRC0;
  real<lower = 0> TVGAMMA;
  real<lower = 0> TVALPHA;
  
  real<lower = 0> omega_mtt;
  real<lower = 0> omega_circ0;
  real<lower = 0> omega_gamma;
  real<lower = 0> omega_alpha;
  
  corr_matrix[3] R;     // Correlation matrix before transforming to Omega.
                        // Can in theory change this to having inputs for
                        // cor_cl_vc, cor_cl_ka, ... and then construct the 
                        // correlation matrix in transformed data, but it's easy
                        // enough to do in R
                        
  corr_matrix[4] R_pd;  // Correlation matrix before transforming to Omega.
                        // Can in theory change this to having inputs for
                        // cor_mtt_circ0, cor_mtt_gamma, ... and then construct 
                        // the correlation matrix in transformed data, but it's 
                        // easy enough to do in R
  
  real<lower = 0> sigma_p;
  real<lower = 0> sigma_a;
  real<lower = -1, upper = 1> cor_p_a;
  
  real<lower = 0> sigma_pd;
  
  real<lower = 0> t_1;
  real<lower = t_1> t_2;
  
}
transformed data{
  
  int n_random = 3;
  int n_random_pd = 4;
  
  int n_cmt = 2;       // number of ODEs in PK model (depot, central)
  int n_cmt_pd = 5;    // number of ODEs in PD system
  int n_cmt_extra = 5; // number of ODEs for AUC, Tmax, Rmin, ...

  vector[n_random] omega = [omega_cl, omega_vc, omega_ka]';
                               
  vector[n_random_pd] omega_pd = [omega_mtt, omega_circ0, omega_gamma, 
                                  omega_alpha]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);
  
  vector[2] sigma = [sigma_p, sigma_a]';
  matrix[2, 2] R_Sigma = rep_matrix(1, 2, 2);
  R_Sigma[1, 2] = cor_p_a;
  R_Sigma[2, 1] = cor_p_a;
  
  matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
  
  matrix[n_random_pd, n_random_pd] L_pd = cholesky_decompose(R_pd);
  
  array[n_cmt + n_cmt_pd + n_cmt_extra] real bioav = 
                                 rep_array(1.0, n_cmt + n_cmt_pd + n_cmt_extra); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt + n_cmt_pd + n_cmt_extra] real tlag = 
                                 rep_array(0.0, n_cmt + n_cmt_pd + n_cmt_extra); // Hardcoding, but could be data or a parameter in another situation
  
  array[1, 2] real x_r = {{t_1, t_2}};
  
}
model{
  
}
generated quantities{
  
  vector[n_total] ipred; // concentration or response with no residual error
  vector[n_total] dv;    // concentration or response with residual error
  
  vector[n_total] auc;            // AUC 
  vector[n_subjects] auc_t1_t2;   // AUC from t1 up to t2
  vector[n_subjects] c_max;       // Cmax
  vector[n_subjects] t_max;       // Tmax foor PK
  vector[n_subjects] t_half;      // half-life
  vector[n_subjects] r_min;       // Minimum response
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  
  vector[n_subjects] MTT;
  vector[n_subjects] CIRC0;
  vector[n_subjects] GAMMA;
  vector[n_subjects] ALPHA;
  
  {
  
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA});
    
    matrix[n_subjects, n_random] eta;   
    matrix[n_subjects, n_random] theta; 
    
    row_vector[n_random_pd] typical_values_pd = to_row_vector({TVMTT, TVCIRC0, 
                                                               TVGAMMA, TVALPHA});
    
    matrix[n_subjects, n_random_pd] eta_pd;   
    matrix[n_subjects, n_random_pd] theta_pd; 
  
    matrix[n_total, n_cmt + n_cmt_pd + n_cmt_extra] x_ipred;
    
    for(i in 1:n_subjects){
      eta[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                           diag_pre_multiply(omega, L))';
                                           
      eta_pd[i, ] = multi_normal_cholesky_rng(rep_vector(0, n_random_pd),
                                              diag_pre_multiply(omega_pd, L_pd))';                                           
    }
    theta = (rep_matrix(typical_values, n_subjects) .* exp(eta));
    theta_pd = (rep_matrix(typical_values_pd, n_subjects) .* exp(eta_pd));

    CL = col(theta, 1);
    VC = col(theta, 2);
    KA = col(theta, 3);
    
    MTT = col(theta_pd, 1);
    CIRC0 = col(theta_pd, 2);
    GAMMA = col(theta_pd, 3);
    ALPHA = col(theta_pd, 4);
    
    for(j in 1:n_subjects){

      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_1cmt_friberg_ode,
                       n_cmt + n_cmt_pd + n_cmt_extra,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {CL[j], VC[j], KA[j], 
                          MTT[j], CIRC0[j], GAMMA[j], ALPHA[j]},
                       bioav, tlag, x_r)';

      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){
          ipred[k] = x_ipred[k, 2] / VC[j];
        }else if(cmt[k] == 4){
          ipred[k] = x_ipred[k, 7] + CIRC0[j];
        }
        auc[k] = x_ipred[k, 8] / VC[j];
      }

      auc_t1_t2[j] = max(x_ipred[subj_start[j]:subj_end[j], 9]) / VC[j];
      c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 10]);
      t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 11]) - t_1;
      t_half[j] = log(2)/(CL[j]/VC[j]);
      
      r_min[j] = min(x_ipred[subj_start[j]:subj_end[j], 12]) + CIRC0[j];
      
    }
    
    for(i in 1:n_total){
      if(ipred[i] == 0){
         dv[i] = 0;
      }else{
        
        if(cmt[i] == 2 || cmt[i] == 4){
          real ipred_tmp = ipred[i];
          real sigma_tmp = cmt[i] == 2 ? 
            sqrt(square(ipred_tmp) * Sigma[1, 1] + Sigma[2, 2] + 2*ipred_tmp*Sigma[2, 1]) :
            sigma_pd;
          
          dv[i] = cmt[i] == 2 ? normal_lb_rng(ipred_tmp, sigma_tmp, 0.0) :
            lognormal_rng(log(ipred_tmp), sigma_tmp);
          
        }
      }
    }
  }
}


