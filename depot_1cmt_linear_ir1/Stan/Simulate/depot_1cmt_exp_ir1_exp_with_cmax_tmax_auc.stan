// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model with Indirect Response 1 PD Model
// IIV on CL, VC, and KA (full covariance matrix)
// IIV on KIN, KOUT, IC50 (full covariance matrix) for PD. IMAX and HILL are 
//   fixed to be 1
// exponential error for PK - DV = IPRED*exp(eps)
// exponential error for PD - DV = IPRED*exp(eps_pd)
// Users choice of general or coupled ODE solution using Torsten
// PK output includes individual Cmax over the whole time period, Tmax between 
//   t1 and t2, AUC since 0 for every timepoint, and AUC between t1 and t2 (like 
//   a dosing interval)
// PD output includes individual Rmin (response decreases to a min), and Tmin
//   betweem t1 and t2 (like the time within a dosing interval at which the 
//   minimum is reached)

functions{
  
  vector depot_1cmt_ir1_ode(real t, vector y, array[] real params, 
                            array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real kin = params[4];
    real kout = params[5];
    real ic50 = params[6];
    real imax = params[7];   // It's fixed to 1 in this particular model
    real hill = params[8];   // It's fixed to 1 in this particular model
    real r_0 = params[9];
    
    real t_1 = x_r[1];
    real t_2 = x_r[2];
    
    real ke = cl/vc;
    
    real conc = y[2]/vc;
    
    real inh = (imax*pow(conc, hill))/(pow(ic50, hill) + pow(conc, hill));
    real response = y[3] + r_0;
    
    real slope_pk = ka*y[1] - ke*y[2];
    real x_pk = slope_pk > 0 && conc > y[6] ? slope_pk/vc : 0;
    real z_pk = t <= t_1 || (slope_pk > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
    
    real slope_pd = kin*(1 - inh) - kout*response;
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
    dydt[8] = x_pd;                                // R_min
    dydt[9] = z_pd;                                // t_min for PD
    
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
  real<lower = 0> TVIC50;
  
  real<lower = 0> omega_kin;
  real<lower = 0> omega_kout;
  real<lower = 0> omega_ic50;
  
  corr_matrix[3] R;     // Correlation matrix before transforming to Omega.
                        // Can in theory change this to having inputs for
                        // cor_cl_vc, cor_cl_ka, ... and then construct the 
                        // correlation matrix in transformed data, but it's easy
                        // enough to do in R
                        
  corr_matrix[3] R_pd;  // Correlation matrix before transforming to Omega.
                        // Can in theory change this to having inputs for
                        // cor_kin_kout, cor_kin_ic50, ... and then construct the 
                        // correlation matrix in transformed data, but it's easy
                        // enough to do in R
  
  real<lower = 0> sigma;
  real<lower = 0> sigma_pd;
  
  real<lower = 0> t_1;
  real<lower = t_1> t_2;
  
}
transformed data{
  
  int n_random = 3;
  int n_random_pd = 3;
  
  int n_cmt = 2;       // number of ODEs in PK model (depot, central)
  int n_cmt_pd = 1;    // number of ODEs in PD system
  int n_cmt_extra = 6; // number of ODEs for AUC, Tmax, Rmin, ...

  vector[n_random] omega = [omega_cl, omega_vc, omega_ka]';
                               
  vector[n_random_pd] omega_pd = [omega_kin, omega_kout, omega_ic50]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);
  
  matrix[n_random_pd, n_random_pd] L_pd = cholesky_decompose(R_pd);
  
  array[n_cmt + n_cmt_pd + n_cmt_extra] real bioav = 
                                 rep_array(1.0, n_cmt + n_cmt_pd + n_cmt_extra); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt + n_cmt_pd + n_cmt_extra] real tlag = 
                                 rep_array(0.0, n_cmt + n_cmt_pd + n_cmt_extra);  // Hardcoding, but could be data or a parameter in another situation
  
  real imax = 1.0;
  real hill = 1.0;
  
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
  vector[n_subjects] t_min;       // Tmin for PD
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  
  vector[n_subjects] KIN;
  vector[n_subjects] KOUT;
  vector[n_subjects] IC50;
  
  {
  
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA});
    
    matrix[n_subjects, n_random] eta;   
    matrix[n_subjects, n_random] theta; 
    
    row_vector[n_random_pd] typical_values_pd = to_row_vector({TVKIN, TVKOUT, 
                                                               TVIC50});
    
    matrix[n_subjects, n_random_pd] eta_pd;   
    matrix[n_subjects, n_random_pd] theta_pd; 
  
    matrix[n_total, n_cmt + n_cmt_pd + n_cmt_extra] x_ipred;
    
    vector[n_subjects] r_0;
    
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
    
    KIN = col(theta_pd, 1);
    KOUT = col(theta_pd, 2);
    IC50 = col(theta_pd, 3);
    
    r_0 = KIN ./ KOUT; 
    
    for(j in 1:n_subjects){

      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(depot_1cmt_ir1_ode,
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
                        KIN[j], KOUT[j], IC50[j], imax, hill, r_0[j]},
                       bioav, tlag, x_r)';

      for(k in subj_start[j]:subj_end[j]){
        if(cmt[k] == 2){
          ipred[k] = x_ipred[k, 2] / VC[j];
        }else if(cmt[k] == 3){
          ipred[k] = x_ipred[k, 3] + r_0[j];
        }
        auc[k] = x_ipred[k, 4] / VC[j];
      }
      
      auc_t1_t2[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]) / VC[j];
      c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 6]);
      t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 7]) - t_1;
      t_half[j] = log(2)/(CL[j]/VC[j]);
      
      r_min[j] = min(x_ipred[subj_start[j]:subj_end[j], 8]) + r_0[j];
      t_min[j] = max(x_ipred[subj_start[j]:subj_end[j], 9]) - t_1;
      
    }
    

    for(i in 1:n_total){
      if(ipred[i] == 0){
         dv[i] = 0;
      }else{
        
        if(cmt[i] == 2 || cmt[i] == 3){
          real log_ipred_tmp = log(ipred[i]);
          real sigma_tmp = cmt[i] == 2 ? sigma : sigma_pd;
          
          dv[i] = lognormal_rng(log_ipred_tmp, sigma_tmp);
          
        }
      }
    }
  }
}


