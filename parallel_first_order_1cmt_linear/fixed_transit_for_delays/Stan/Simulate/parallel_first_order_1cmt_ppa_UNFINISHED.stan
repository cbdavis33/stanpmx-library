// One-compartment PK Model with Parallel absorption
// This is flexible so that any number of parallel absorption processes can take 
//   place (n_depot >= 2).
// A delay in absorption for each process is implemented with transit 
//   compartments (n_transit). This number of transit compartments must be 
//   constant across depots, with the exception being that the fastest process
//   can have no delay if desired. Control this with "delay_on_fastest_depot".
// There will be n_depot KA values for each subject. i.e., TVKA and KAi are 
//   vectors. KA is intended to be ordered so that 
//   KA_1 < KA_2 < ... < KA_n_depot. I won't do this exactly, but I will have
//   TVKA_1 < TVKA_2 < ... < TVKA_n_depot and hope the individual effects don't
//   mess that up
// There will be n_depot MTT values for each subject. i.e., TVMTT and MTTi are 
//   vectors. MTT is intended to be ordered so that 
//   MTT_1 > MTT_2 > ... > MTT_n_depot. I won't do this exactly, but I will have
//   TVMTT_1 > TVMTT_2 > ... > TVMTT_n_depot and hope the individual effects
//   don't mess that up
// FRAC is a simplex of length n_depot that tells how much of the dose goes into
//   each absorption process
// IIV on CL, VC, KA, MTT (full covariance matrix)

// Either of matrix-exponential or general ODE solution using Torsten
// proportional plus additive error - DV = IPRED*(1 + eps_p) + eps_a
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   generates values from a normal that is truncated below at 0

functions{
  
  real normal_lb_rng(real mu, real sigma, real lb){
     
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y;

  }
  
  // vector transit_fixed_1cmt_ode(real t, vector y, array[] real params, 
  //                               array[] real x_r, array[] int x_i){
  //   
  //   real cl = params[1];
  //   real vc = params[2];
  //   real ka = params[3];
  //   real mtt = params[4];
  //   
  //   int n_transit = x_i[1];
  //   
  //   real ktr = (n_transit + 1)/mtt;
  //   
  //   real ke = cl/vc;
  //   
  //   vector[n_transit + 3] dydt;
  // 
  //   dydt[1] = -ktr*y[1];                               // 0th transit compartment (depot)
  //   for(i in 2:(n_transit + 1)){
  //     dydt[i] = ktr*(y[i - 1] - y[i]);
  //   }
  //   dydt[n_transit + 2] = ktr*y[n_transit + 1] - ka*y[n_transit + 2];
  //   dydt[n_transit + 3] = ka*y[n_transit + 2] - ke*y[n_transit + 3];
  //   
  //   return dydt;
  // }

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
  
  int<lower = 2> n_depots;
  int<lower = 0, upper = 1> delay_on_fastest_depot;
  int<lower = 1> n_transit;
  
  real<lower = 0> TVCL;
  real<lower = 0> TVVC;
  positive_ordered[n_depots] TVKA;
  array[n_depots] real TVMTT; // SHOULD BE DECREASING ORDERED
  
  simplex[n_depots] TVFRAC;
  
  real<lower = 0> omega_cl;
  real<lower = 0> omega_vc;
  vector<lower = 0>[n_depots] omega_ka;
  vector<lower = 0>[n_depots] omega_mtt;
  
  corr_matrix[n_depots*2 + 2] R;  // Correlation matrix before transforming to Omega.
                                  // Can in theory change this to having inputs for
                                  // cor_cl_vc, cor_cl_ka, ... and then construct the 
                                  // correlation matrix in transformed data, but it's easy
                                  // enough to do in R
  
  real<lower = 0> sigma_p;
  real<lower = 0> sigma_a;
  real<lower = -1, upper = 1> cor_p_a;
  
  int<lower = 1, upper = 2> solver; // 1 = matrix-exponential, 2 = rk45
  
}
transformed data{
  
  int n_random = n_depots*2 + 2;           // Number of random effects

  int n_cmt;
  if(delay_on_fastest_depot == 1){
    n_cmt = n_depots*(n_transit + 2) + 1;
    n_depots_with_delay = n_depots;
  }else{
    n_cmt = (n_depot - 1)*(n_transit + 2) + 2;
    n_depots_with_delay = n_depots - 1;
    for(i in 1:n_total){
      if(cmt[i] == n_depot){
        cmt[i] = n_cmt - 1; // TODO: Might not be able to overwrite cmt like this
      }
    }
  }

  // TODO: Figure out this omega = ... line
  vector[n_random] omega = append_row(append_row(to_vector(omega_ka), omega_cl), omega_vc);
  
  matrix[n_random, n_random] L = cholesky_decompose(R);

  vector[2] sigma = [sigma_p, sigma_a]';
  matrix[2, 2] R_Sigma = rep_matrix(1, 2, 2);
  R_Sigma[1, 2] = cor_p_a;
  R_Sigma[2, 1] = cor_p_a;
  
  matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);  // Hardcoding, but could be data or a parameter in another situation
  
  array[1, 2] real x_r = {{8675309, 5555555}}; // This is a placeholder of nonsense so I can put it in the Torsten function
  
  // TODO: figure out which of the is correct
  array[1, 1] int x_i = {{n_depots_with_delay, n_transit}};
  array[1, 1] int x_i = {{n_depots, n_transit}};
  
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
        K[(n_transit + 3), (n_transit + 3)] = -(CL[j]/VC[j]);

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
        real ipred_tmp = ipred[i];
        real sigma_tmp = sqrt(square(ipred_tmp) * Sigma[1, 1] + Sigma[2, 2] + 
                              2*ipred_tmp*Sigma[2, 1]);
        dv[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
        
      }
    }
  }
}

