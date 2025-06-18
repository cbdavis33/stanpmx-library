// First Order Absorption (oral/subcutaneous)
// One-compartment PK Model
// IIV on CL, VC, and Ka (full covariance matrix)
// proportional plus additive error with a Student's-t distribution - 
//   DV = IPRED*(1 + eps_p) + eps_a
// User's choice of analytical solution or general ODE solution
// General ODE solution using Torsten will get out individual estimates of AUC, 
//   Cmax, Tmax, ... Analytical will not
// Predictions are generated from a Student's-t that is truncated below at 0

functions{
  
  vector n_t_q_inner_coeffs(int n, real nu) {
    if (n <= 0)
      reject("n must be positive; found n = ", n);
  
  real g = sqrt(nu/2) * exp(lgamma(nu/2) - lgamma((nu + 1)/2));
    vector[3] x0 = [
      g,
      ((nu + 1)*g^3 - nu*g)/(6*nu),
      ((7*nu^2 + 8*nu + 1)*g^5 - 10*(nu^2 + nu)*g^3 + 3*nu^2*g)/(120*nu^2)
    ]';
  
    if(n < 4)
      return x0[1:n];
  
    vector[n] x = append_row(x0, rep_vector(0, n - 3));
  
    for(i in 2:(n - 2)) {
      x[i + 2] = -(2*i + 1)*x[i + 1];
      
      for(l in 0:i) {
        for(m in 0:(i - l)) {
          x[i + 2] = x[i + 2] + ((1 + 1/nu)*(2*l + 1)*(2*m + 1) - 
                       2*m*(2*m + 1)/nu)*x[i - l - m + 1]*x[l + 1]*x[m + 1];
        }
      }
      for (l in 0:(i - 1)) {
        for (m in 0:(i - l - 1)) {
          x[i + 2] = x[i + 2] - (1/nu)*(2*m + 1)*x[i - l - m]*x[l + 1]*x[m + 1];
        }
      }
      
      x[i + 2] = x[i + 2]/(2*i + 3)/(2*i + 2);
    }
  
    return x;
}

  real n_t_q(real z, real nu) {
    if(abs(z) < 2*sqrt(nu)) {
      vector[11] a = n_t_q_inner_coeffs(11, nu);
      for (i in 0:10)
        a[i + 1] *= z^(2*i + 1);
      return sum(a);
    } else {
      real w = (1 - std_normal_cdf(abs(z)))*nu*sqrt(pi())*exp(lgamma(nu/2) - lgamma((nu + 1)/2));
      real sign_z = z > 0 ? 1 : z < 0 ? -1 : 0;
      return sign_z*sqrt(nu)*w^(-1/nu)*(1 - (nu + 1)/(2*(nu + 2))*w^(2/nu));
    }
  }

  real student_t_alternate_rng(real nu, real mu, real sigma) {
    return mu + sigma * n_t_q(std_normal_rng(), nu);
  }

  real student_t_lub_rng(real nu, real mu, real sigma, real lb, real ub) {
    real p_lb = student_t_cdf(lb | nu, mu, sigma);
    real p_ub = student_t_cdf(ub | nu, mu, sigma);
    real u = uniform_rng(p_lb, p_ub);
    real z = inv_Phi(u);
    return mu + sigma * n_t_q(z, nu);
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
    
    vector[6] dydt;

    dydt[1] = -ka*y[1];                            // depot
    dydt[2] = slope;                               // central
    dydt[3] = y[2];                                // AUC
    dydt[4] = t >= t_1 && t <= t_2 ? y[2] : 0;     // AUC_t_1-t_2
    dydt[5] = x;                                   // C_max
    dydt[6] = z;                                   // t_max
    
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
  
  int n_random = 3; // Number of random effects
  int n_cmt = want_auc_cmax ? 6 : 2; // Number of compartments - depot, central (AUC, AUC_ss, Cmax_ss, Tmax_ss))
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt);
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
  array[1, 2] real x_r = {{t_1, t_2}};

}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = TVCL/TVVC> TVKA;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  vector<lower = 0>[2] sigma;
  cholesky_factor_corr[2] L_Sigma;
  
  real nu;
  
  matrix[n_random, n_subjects] Z;
  
}
generated quantities{
  
  vector[n_time_new] pred;       // f(TVs, x, eta = 0) 
  vector[n_time_new] epred_stan; // f(TVs, x, eta = eta_new), eta_new ~ multi_normal(0, Omega) 
  vector[n_time_new] epred;      // epred_stan + error
  vector[n_time_new] ipred;      // f(TVs, x, eta = eta_i), eta_i are etas for observed subjects
  vector[n_time_new] dv;         // ipred + error
  
  vector[want_auc_cmax ? n_time_new : 0] auc;     // AUC from 0 up until t
  vector[want_auc_cmax ? n_subjects : 0] auc_ss;  // AUC from t1 up to t2 (AUC_ss)
  vector[want_auc_cmax ? n_subjects : 0] c_max;   // Cmax between t1 and t2 (c_max_ss)
  vector[want_auc_cmax ? n_subjects : 0] t_max;   // Tmax between t1 and t2, then subtract off t1
  vector[n_subjects] t_half;                      // half-life
  
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] KA;
  vector[n_subjects] KE;
 
  {
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVKA});
    
    matrix[2, 2] R_Sigma = multiply_lower_tri_self_transpose(L_Sigma);
    matrix[2, 2] Sigma = quad_form_diag(R_Sigma, sigma);
    
    real sigma_sq_p = Sigma[1, 1];
    real sigma_sq_a = Sigma[2, 2];
    real sigma_p_a = Sigma[1, 2];

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

    matrix[n_time_new, 2] x_pred;
    matrix[n_time_new, 2] x_epred;
    matrix[n_time_new, n_cmt] x_ipred;
    
    for(j in 1:n_subjects){
      
      row_vector[n_random] theta_j = theta[j]; // access the parameters for subject j
      row_vector[n_random] theta_j_new = theta_new[j]; // access the parameters for subject j's epred
      
      real cl_p = TVCL;
      real vc_p = TVVC;
      real ka_p = TVKA;
      
      CL[j] = theta_j[1];
      VC[j] = theta_j[2];
      KA[j] = theta_j[3];
      KE[j] = CL[j]/VC[j];
      t_half[j] = log(2)/KE[j];
      
      CL_new[j] = theta_j_new[1];
      VC_new[j] = theta_j_new[2];
      KA_new[j] = theta_j_new[3];
      
      // x_pred[subj_start[j]:subj_end[j],] =
      //   pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
      //                    amt[subj_start[j]:subj_end[j]],
      //                    rate[subj_start[j]:subj_end[j]],
      //                    ii[subj_start[j]:subj_end[j]],
      //                    evid[subj_start[j]:subj_end[j]],
      //                    cmt[subj_start[j]:subj_end[j]],
      //                    addl[subj_start[j]:subj_end[j]],
      //                    ss[subj_start[j]:subj_end[j]],
      //                    {cl_p, vc_p, ka_p})';
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
      //                    {CL_new[j], VC_new[j], KA_new[j]})';
      
      matrix[2, 2] K_epred = rep_matrix(0, 2, 2);
      matrix[2, 2] K_p = rep_matrix(0, 2, 2);
      
      K_epred[1, 1] = -KA_new[j];
      K_epred[2, 1] = KA_new[j];
      K_epred[2, 2] = -CL_new[j]/VC_new[j];
      
      x_epred[subj_start[j]:subj_end[j], ] =
        pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         K_epred, bioav, tlag)';
                         
      epred_stan[subj_start[j]:subj_end[j]] =
        x_epred[subj_start[j]:subj_end[j], 2] ./ VC_new[j];
                         
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
                         K_p, bioav, tlag)';
                         
      pred[subj_start[j]:subj_end[j]] =
        x_pred[subj_start[j]:subj_end[j], 2] ./ vc_p;
      
      
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
                         {CL[j], VC[j], KA[j]}, bioav, tlag, x_r)';
                        
        ipred[subj_start[j]:subj_end[j]] = 
          x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
        
        auc[subj_start[j]:subj_end[j]] = 
                                x_ipred[subj_start[j]:subj_end[j], 3] ./ VC[j];
        auc_ss[j] = max(x_ipred[subj_start[j]:subj_end[j], 4]) / VC[j];
        c_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 5]);
        t_max[j] = max(x_ipred[subj_start[j]:subj_end[j], 6]) - t_1;
          
      }else{
        
        // x_ipred[subj_start[j]:subj_end[j],] =
        //   pmx_solve_onecpt(time[subj_start[j]:subj_end[j]],
        //                    amt[subj_start[j]:subj_end[j]],
        //                    rate[subj_start[j]:subj_end[j]],
        //                    ii[subj_start[j]:subj_end[j]],
        //                    evid[subj_start[j]:subj_end[j]],
        //                    cmt[subj_start[j]:subj_end[j]],
        //                    addl[subj_start[j]:subj_end[j]],
        //                    ss[subj_start[j]:subj_end[j]],
        //                    {CL[j], VC[j], KA[j]})';
                           
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
                           K, bioav, tlag)';

        ipred[subj_start[j]:subj_end[j]] = 
          x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
        
      }

    }
  
    for(i in 1:n_time_new){
      
      if(epred_stan[i] == 0){
        epred[i] = 0;
      }else{
        real epred_tmp = epred_stan[i];
        real sigma_tmp_e = sqrt(square(epred_tmp) * sigma_sq_p + sigma_sq_a + 
                                2*epred_tmp*sigma_p_a);
        epred[i] = student_t_lub_rng(nu, epred_tmp, sigma_tmp_e, 0.0, 
                                     positive_infinity());
      }
      
      if(ipred[i] == 0){
        dv[i] = 0;
      }else{
        real ipred_tmp = ipred[i];
        real sigma_tmp = sqrt(square(ipred_tmp) * sigma_sq_p + sigma_sq_a + 
                              2*ipred_tmp*sigma_p_a);
        dv[i] = student_t_lub_rng(nu, ipred_tmp, sigma_tmp, 0.0, 
                                  positive_infinity());
      } 
      
    }
  }
}
