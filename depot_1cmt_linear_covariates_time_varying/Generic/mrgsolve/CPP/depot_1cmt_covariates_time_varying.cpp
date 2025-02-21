$PARAM @annotated
TVCL  : 0.5 : Population clearance  
TVVC  :  12 : Population Central Volume
TVKA  :  1  : Population absorption 
  
THETA_CL_WT    : 0.75 : body weight effect on clearance
THETA_VC_WT    : 1    : body weight effect on central volume
THETA_KA_CMPPI :  0   : cmppi effect on absorption
THETA_CL_EGFR  :  0   : eGFR effect on clearance
  
WT    :  70 : Body weight
EGFR  :  90 : eGFR
CMPPI :  0  : protein pump inhibitor
  
  
$SET delta=0.1, end=240

$CMT @annotated
DEPOT : Depot
CENT  : Central compartment

$MAIN
double CLi = exp(log(TVCL) + THETA_CL_WT*log(WT/70) + THETA_CL_EGFR*log(EGFR/90) + ECL);
double VCi = exp(log(TVVC) + THETA_VC_WT*log(WT/70) + EVC);
double KAi = exp(log(TVKA) + THETA_KA_CMPPI*CMPPI + EKA);

$OMEGA @name IIV @labels ECL EVC EKA
0 0 0

$SIGMA 0
  
$ODE
dxdt_DEPOT = -KAi*DEPOT;
dxdt_CENT = KAi*DEPOT - (CLi/VCi)*CENT;

$CAPTURE CLi VCi KAi
  
$TABLE
capture IPRED = CENT/VCi;
