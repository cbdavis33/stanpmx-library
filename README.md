# A Library of Stan Models and R Scripts for PMx

Note: The workflow is a little out-of-date for most of these repositories. As of June 7, 2025, the following are fully updated

-   depot_1cmt_linear

-   depot_1cmt_linear_covariates

-   depot_2cmt_linear

-   iv_1cmt_linear

-   iv_2cmt_linear

-   iv_2cmt_linear_covariates

-   bioav_iv_and_oral_1cmt_linear

-   bioav_iv_and_oral_2cmt_linear

The others have correct Stan code, so if you just want to fit, they're still good, but I haven't updated the post-processing and analyses you might want to do after fitting, yet. I'll be going through them and updating them as time allows.

This is a collection of pharmacometric models written in Stan + Torsten and R. In this repo are templates for many common and less-common models in Stan + Torsten - basic one- and two-compartment IV and/or oral models, Michaelis-Menten elimination, indirect response models, effect-compartment models, the Friberg-Karlsson neutropenia model, the Savic transit compartment absorption model, delayed absorption models with a fixed number of transit compartments, models with covariates (some time-varying), BLOQ-handling, exposure-response analyses (eventually), etc. …. [bayespmx.github.io](https://bayespmx.github.io) (incomplete and in the early stages, but updated periodically) and [stanpmx.github.io](https://stanpmx.github.io) (no longer maintained) provide some tutorials and guidance. I believe that making these models public along with some tutorials will give users in the PMx community enough of a knowledge platform and a code template that they can implement a fully Bayesian workflow that involves fake data simulation, model fitting, post-processing of the posterior distribution, model diagnostics, model selection, and making predictions and future simulations, for any model that they desire.

Hopefully, these resources will help aspiring Bayesians overcome some of the early barriers to entry (knowledge, code) while providing a collection of scripts and models for the PMx community to use as a stepping stone to implement their own Bayesian models in their drug-development process. It will hopefully become a user-driven repository where pharmacometricians can submit their own models and tutorials that they think others in the PMx world could benefit from.

Here's the rough idea - each directory I hope is named so that it makes it obvious what model is in there, e.g. iv_1cmt_linear is for a one-compartment model with IV infusion and linear elimination, transit_savic_2cmt_linear is a two-compartment model with an absorption delay modeled by the [Savic transit compartment model](https://pubmed.ncbi.nlm.nih.gov/17653836/). Each directory should have models with lognormal (\*\_exp), proportional (\*\_prop), and proportional-plus-additive error (\*\_ppa). I hope it's easy to understand what each file does by the combination of file path and filename. e.g. iv_1cmt_linear/Stan/Fit/iv_1cmt_prop.stan is the Stan model that fits one-cmt IV data with a proportional error model, and iv_1cmt_linear/R/Fit/iv_1cmt_prop.R is the R code needed to import the NONMEM dataset, wrangle it into a Stan-ready format, sample from the posterior, and save the fitted object.

Basically, each directory is meant to be able to

1)  \*/\*/Simulate/ - Simulate fake data. Roughly a substitute for mrgsolve or RxODE. It's not necessarily as flexible and comprehensive as those, but it's a good part of the workflow, partly so you can create some fake data to fit, and partly so you can quickly simulate data to check you've written your ODEs correctly. I think it's an important step in the process. Since proportional error is a special case of proportional-plus-additive error (set sigma_a and cor_p_a to 0), there's only \*\_exp and \*\_ppa. Many of the \*/Simulate directories have something like \*\_with_auc.\*. These simulate the data with a system of ODEs plus an AUC compartment in case that's a desired output. Some of them have a Cmax compartment, too. Those models that have an analytical solution and/or are linear ODEs have files without \_with_auc. These simulation files give the option of simulating using the analytical, matrix-exponential, or ODE solution.

2)  \*/\*/Fit/ - Fit the model and save the fitted object. Most of the models implement [within-chain parallelization with reduce_sum()](https://mc-stan.org/docs/stan-users-guide/reduce-sum.html) with [some code](https://bayespmx.github.io/tutorials/Threading-for-Within-Chain-Parallelization.html#example-one-compartment-iv) I wrote. The ones that don't do this have \*\_no_threading in the filename. The no-threading files are mostly for speed testing. I wouldn't bother actually fitting with those unless you're on a very basic machine with 8 cores or fewer and the model is extremely simple. The threaded models are faster than the non-threaded models if you have a computing infrastructure (high-end laptop, desktop. or HPC) that will allow for more than 2 threads per chain. In theory, there are cases where threading isn't faster and is actually slower, but not for any PopPK/PopPKPD models. Many of the models also have a \*\_prop_all_solvers.\* file. This one is meant to fit the proportional error model while giving the user the option to use the analytical, matrix-exponential, or ODE solution. I always wrote this one to test the speed of the solvers. For some models, the analytical solution is fastest (e.g. depot_2cmt_linear), and for some models, the matrix-exponential solution was fastest (e.g. iv_2cmt_linear). From there, I wrote the rest of the models with whatever was fastest. Also, I mostly use traditional NONMEM/PMx names (like TVCL vs. the CLhat that Metrum uses). Also, these files spit out a lot of quantities people like seeing, e.g. IWRES, PRED/EPRED/IPRED, ..., and it names the outputs something that make sense, so I have something like omega_cl, omega_vc, omega_q, omega_vp rather than omega[1], omega[2], omega[3], omega[4] and cor_cl_vc rather than R[1, 2]... Hopefully this all makes sense, and you can figure out what everything is by its name in the code.

3)  \*/Stan/Fits/ - Save the fitted objects here (from \*/R/Fit/\*.R). There is also an Output/ subdirectory here where I write the CSVs so that I can check on the progress for long-running models. There is also a Stan_Data/ subdirectory here where I write the stan_data list that goes into the fit.

4)  \*/R/Post_Process/ - There are files for looking at posterior summaries, MCMC diagnostics, and model diagnostics like, DV vs. PRED/EPRED/IPRED, individual plots, VPCs, NPDEs, residuals, shrinkage, .... \*\_compare_prior_and_posterior.\* also reads in the fitted object, simulates from the prior, then compares the prior and the posterior.

5)  \*/\*/Predict/ - \*\_predict_observed_subjects.\* simulates from the posterior for the observed subjects. Basically meaning it takes the actual observed subjects, their actual dosing regimen, and simulates from their individual posterior at the timepoints you want to see. This helps make prettier PPC pictures if you tell it to simulate at a denser grid than the observations. Ig all you want to do is make nicer pictures, then you can do these prediction with whatever is the fastest solver. If you want other outputs about the observed subjects like AUC, Cmax, Tmax, and AUC-between-time1-and-time2, then you can use the ODE solution - it's written with extra compartments to get model-based estimates for these quantities. \*\_predict_new_subjects.\* simulates from the posterior for new subjects. You can choose the dosing regimen (and covariates if applicable) and then simulate new subjects. The difference here from \*\_predict_observed_subjects.\* is the simulation of new etas for the new subjects. Pretty standard. The main use for this is simulating potential new dosing regimens.

6)  \*/R/VPC/ - Do a VPC and PCVPC using the vpc package.

7)  \*/R/NPDE/ - NPDEs using the npde package.

8)  \*/<Stan or R>/Cross_Validation/ - Cross-validation. Leave-one-out (LOO) cross-validation is done with the fitting. The .stan files here perform Leave-one-group-out (LOGO) cross-validation
