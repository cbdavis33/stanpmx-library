// One-compartment PK Model with Transit Compartment absorption
// IIV on CL, VC, KA, MTT (full covariance matrix)
// n_transit is fixed to a positive integer that can be defined as data
// The 0th transit compartment is where the dosing happens (cmt = 1). This could
//   also be called the Depot
// Either of matrix-exponential or general ODE solution using Torsten
// exponential error - DV = IPRED*exp(eps)

functions{
  
  vector transit_fixed_1cmt_ode(real t, vector y, array[] real params, 
                                array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real ka = params[3];
    real mtt = params[4];
    
    int n_transit = x_i[1];
    
    real ktr = (n_transit + 1)/mtt;
    
    real ke = cl/vc;
    
    vector[n_transit + 3] dydt;

    dydt[1] = -ktr*y[1];                               // 0th transit compartment (depot)
    for(i in 2:(n_transit + 1)){
      dydt[i] = ktr*(y[i - 1] - y[i]);
    }
    dydt[n_transit + 2] = ktr*y[n_transit + 1] - ka*y[n_transit + 2];
    dydt[n_transit + 3] = ka*y[n_transit + 2] - ke*y[n_transit + 3];
    
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
  real<lower = 0> TVMTT;
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  real<lower = 0> omega_ka;
  real<lower = 0> omega_mtt;
  
  corr_matrix[4] R;  // Correlation matrix before transforming to Omega.
                     // Can in theory change this to having inputs for
                     // cor_cl_vc, cor_cl_ka, ... and then construct the 
                     // correlation matrix in transformed data, but it's easy
                     // enough to do in R
  
  real<lower = 0> sigma;

  int<lower = 1> n_transit;
  
  int<lower = 1, upper = 2> solver; // 1 = matrix-exponential, 2 = rk45
  
}
transformed data{
  
  int n_random = 4;           // Number of random effects
  int n_cmt = n_transit + 3;  // Depot, tr_1, ..., tr_n, absorption, central

  vector[n_random] omega = [omega_cl, omega_vc, omega_ka, omega_mtt]';
  
  matrix[n_random, n_random] L = cholesky_decompose(R);
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);  // Hardcoding, but could be data or a parameter in another situation
  
  array[1, 2] real x_r = {{8675309, 5555555}}; // This is a placeholder of nonsense so I can put it in the Torsten function
  
  array[1, 1] int x_i = {{n_transit}};
  
}
model{
  
}
generated quantities{
  
  vector[n_total] ipred;  // concentration with no residual error
  vector[n_total] dv;     // concentration with residual error
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  vector[n_subjects] MTT;
  vector[n_subjects] KTR;
  
  {
    
    vector[n_random] typical_values = to_vector({TVCL, TVVC, TVKA, TVMTT});
    
    matrix[n_random, n_subjects] eta;   
    matrix[n_subjects, n_random] theta; 
  
    matrix[n_total, n_cmt] x_ipred;
    
    for(i in 1:n_subjects){
      eta[, i] = multi_normal_cholesky_rng(rep_vector(0, n_random),
                                           diag_pre_multiply(omega, L));
    }
    theta = (rep_matrix(typical_values, n_subjects) .* exp(eta))';

    CL = col(theta, 1);
    VC = col(theta, 2);
    KA = col(theta, 3);
    MTT = col(theta, 4);
    KTR = (n_transit + 1) ./ MTT;
  
    for(j in 1:n_subjects){
      
      if(solver == 1){
        
        matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
      
        for(i in 1:(n_transit + 1)){
          K[i, i] = -KTR[j];
          K[(i + 1), i] = KTR[j];
        }
        
        K[(n_transit + 2), (n_transit + 2)] = -KA[j];
        K[(n_transit + 3), (n_transit + 2)] = KA[j];
        K[(n_transit + 3), (n_transit + 3)] = -CL[j]/VC[j];

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
      }else{
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(transit_fixed_1cmt_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j], KA[j], MTT[j]}, 
                         bioav, tlag, x_r, x_i)';
                         
      }

      ipred[subj_start[j]:subj_end[j]] =
        x_ipred[subj_start[j]:subj_end[j], (n_transit + 3)] ./ VC[j];
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

