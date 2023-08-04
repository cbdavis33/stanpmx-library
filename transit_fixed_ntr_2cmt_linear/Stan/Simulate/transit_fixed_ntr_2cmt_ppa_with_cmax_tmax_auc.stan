// Two-compartment PK Model with Transit Compartment absorption
// IIV on CL, VC, Q, VP, KA, MTT (full covariance matrix)
// n_transit is fixed to a positive integer that can be defined as data
// The 0th transit compartment is where the dosing happens (cmt = 1). This could
//   also be called the Depot
// General ODE solution using Torsten
// proportional plus additive error - DV = IPRED*(1 + eps_p) + eps_a
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
  
  // vector transit_fixed_2cmt_ode(real t, vector y, array[] real params,
  //                               array[] real x_r, array[] int x_i){
  // 
  //   real cl = params[1];
  //   real vc = params[2];
  //   real q = params[3];
  //   real vp = params[4];
  //   real ka = params[5];
  //   real mtt = params[6];
  // 
  //   real t_1 = x_r[1];
  //   real t_2 = x_r[2];
  // 
  //   int n_transit = x_i[1];
  // 
  //   real ktr = (n_transit + 1)/mtt;
  // 
  //   real ke = cl/vc;
  //   real k_cp = q/vc;
  //   real k_pc = q/vp;
  // 
  //   real slope = ka*y[n_transit + 2] - (ke + k_cp)*y[n_transit + 3] +
  //                         k_pc*y[n_transit + 4];
  //   real x = slope > 0 && y[n_transit + 3]/vc > y[n_transit + 7] ? slope/vc : 0;
  //   real z = t <= t_1 || (slope > 0 && t >= t_1 && t <= t_2) ? 1 : 0;
  // 
  //   vector[n_transit + 8] dydt;
  // 
  //   dydt[1] = -ktr*y[1];                               // 0th transit compartment (depot)
  //   for(i in 2:(n_transit + 1)){
  //     dydt[i] = ktr*(y[i - 1] - y[i]);
  //   }
  //   dydt[n_transit + 2] = ktr*y[n_transit + 1] - ka*y[n_transit + 2];
  //   dydt[n_transit + 3] = slope;
  //   dydt[n_transit + 4] = k_cp*y[n_transit + 3] - k_pc*y[n_transit + 4];
  //   dydt[n_transit + 5] = y[n_transit + 3];                    // AUC
  //   dydt[n_transit + 6] = t >= t_1 && t <= t_2 ? y[2] : 0;     // AUC_t_1-t_2
  //   dydt[n_transit + 7] = x;                                   // C_max
  //   dydt[n_transit + 8] = z;                                   // t_max
  // 
  //   return dydt;
  // }
  
    vector transit_fixed_2cmt_ode(real t, vector y, array[] real params,
                                array[] real x_r, array[] int x_i){

    real cl = params[1];
    real vc = params[2];
    real q = params[3];
    real vp = params[4];
    real ka = params[5];
    real mtt = params[6];

    real t_1 = x_r[1];
    real t_2 = x_r[2];

    int n_transit = x_i[1];

    real ktr = (n_transit + 1)/mtt;

    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;

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
    dydt[n_transit + 5] = y[n_transit + 3];                    // AUC
    dydt[n_transit + 6] = t >= t_1 && t <= t_2 ? y[2] : 0;     // AUC_t_1-t_2
    dydt[n_transit + 7] = x;                                   // C_max
    dydt[n_transit + 8] = z;                                   // t_max

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
  real<lower = 0> TVMTT;
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_q;
  real<lower = 0> omega_vp;
  real<lower = 0> omega_ka;
  real<lower = 0> omega_mtt;
  
  corr_matrix[6] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_cl_vc, cor_cl_q, ... and then construct the 
                     // correlation matrix in transformed data, but it's easy
                     // enough to do in R
  
  real<lower = 0> sigma_p;
  real<lower = 0> sigma_a;
  real<lower = -1, upper = 1> cor_p_a;

  int<lower = 1> n_transit;
  
  real<lower = 0> t_1;
  real<lower = 0> t_2;
  
}
transformed data{
  
  int n_random = 6;           // Number of random effects
  int n_cmt = n_transit + 8;  // Depot, tr_1, ..., tr_n, absorption, central, peripheral
                              //   AUC, AUC_t1-t2, Cmax, Tmax

  vector[n_random] omega = [omega_cl, omega_vc, omega_q, omega_vp, omega_ka,
                            omega_mtt]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);

  vector[2] sigma = [sigma_p, sigma_a]';
  matrix[2, 2] R_Sigma = rep_matrix(1, 2, 2);
  R_Sigma[1, 2] = cor_p_a;
  R_Sigma[2, 1] = cor_p_a;
  
  matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);  // Hardcoding, but could be data or a parameter in another situation
  
  array[1, 2] real x_r = {{t_1, t_2}};
  
  array[1, 1] int x_i = {{n_transit}};
  
}
model{
  
}
generated quantities{
  
  vector[n_total] ipred;  // concentration with no residual error
  vector[n_total] dv;     // concentration with residual error
  
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
  vector[n_subjects] MTT;
  vector[n_subjects] KTR;
  
  {
    
    vector[n_random] typical_values = to_vector({TVCL, TVVC, TVQ, TVVP, TVKA, 
                                                 TVMTT});
    
    matrix[n_random, n_subjects] eta;   
    matrix[n_subjects, n_random] theta; 
  
    matrix[n_total, n_cmt] x_ipred;
    
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
    MTT = col(theta, 6);
    KTR = (n_transit + 1) ./ MTT;
    
    alpha = 0.5*(CL./VC + Q./VC + Q./VP + 
                 sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
    beta = 0.5*(CL./VC + Q./VC + Q./VP - 
                 sqrt((CL./VC + Q./VC + Q./VP)^2 - 4*CL./VC.*Q./VP));
  
    for(j in 1:n_subjects){

      x_ipred[subj_start[j]:subj_end[j],] =
        pmx_solve_rk45(transit_fixed_2cmt_ode,
                       n_cmt,
                       time[subj_start[j]:subj_end[j]],
                       amt[subj_start[j]:subj_end[j]],
                       rate[subj_start[j]:subj_end[j]],
                       ii[subj_start[j]:subj_end[j]],
                       evid[subj_start[j]:subj_end[j]],
                       cmt[subj_start[j]:subj_end[j]],
                       addl[subj_start[j]:subj_end[j]],
                       ss[subj_start[j]:subj_end[j]],
                       {CL[j], VC[j], Q[j], VP[j], KA[j], MTT[j]},
                       bioav, tlag, x_r, x_i)';


      ipred[subj_start[j]:subj_end[j]] =
        x_ipred[subj_start[j]:subj_end[j], (n_transit + 3)] ./ VC[j];
        
      auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 4] ./ VC[j];
      
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

