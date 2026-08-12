function output=EZFHorz_CrossTests3_d_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% NEW cross test, specific to Epstein-Zin (no QH analog): a model with markov z plus iid e,
% versus the SAME two shocks set up as a single joint (two-variable) markov with
% n_z=[n_z,n_e] and pi_z=kron(ones(n_e,1)*pi_e',pi_z).
% For vNM/QH this equality is trivially true. For EZ it is the sharpest test that the codes take
% ONE certainty-equivalent over the joint distribution of (zprime,eprime), rather than nesting a
% CE over e inside a CE over z (nested CEs = joint CE only when gamma=1/phi, so this test only
% discriminates because the default ezgamma, ezphi, ezrisk are away from that knife-edge).
% V and StationaryDist agree to summation-order roundoff (~1e-13); Policy agrees exactly
% (barring knife-edge argmax ties).

ReturnFn_ze_cons=@(d,aprime,a,z,e,r,w,kappa_j,varphi,agej,Jr,pension)...
    EZReturnFn_cons_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,varphi,agej,Jr,pension);
ReturnFn_ze_posU=@(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)...
    EZReturnFn_positiveUtils_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
ReturnFn_ze_negU=@(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)...
    EZReturnFn_negativeUtils_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);

% Model B objects: the same (z,e) as a single joint markov (z takes the first slot, e the second,
% so the same ReturnFn works for both models)
n_zB=[n_z,vfoptionsbaseline.n_e];
z_gridB=[z_grid; vfoptionsbaseline.e_grid];
pi_zB=kron(ones(vfoptionsbaseline.n_e,1)*vfoptionsbaseline.pi_e',pi_z);


%% Epstein-Zin case 1: Consumption-units (EZutils=0)
% Model A: markov z plus iid e
vfoptionsA=struct();
vfoptionsA.exoticpreferences='EpsteinZin';
vfoptionsA.EZutils=0;
vfoptionsA.EZriskaversion='ezgamma';
vfoptionsA.EZeis='ezphi';
vfoptionsA.n_e=vfoptionsbaseline.n_e;
vfoptionsA.e_grid=vfoptionsbaseline.e_grid;
vfoptionsA.pi_e=vfoptionsbaseline.pi_e;
simoptionsA=struct();
simoptionsA.n_e=simoptionsbaseline.n_e;
simoptionsA.e_grid=simoptionsbaseline.e_grid;
simoptionsA.pi_e=simoptionsbaseline.pi_e;
jequaloneDistA=zeros(n_a,n_z,vfoptionsA.n_e,'gpuArray');
jequaloneDistA(1,ceil(n_z/2),ceil(vfoptionsA.n_e/2))=1; % no assets, midpoint shocks
[VA,PolicyA]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze_cons,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDistA=StationaryDist_FHorz_Case1(jequaloneDistA,AgeWeightParamNames,PolicyA,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

% Model B: the same two shocks as a single joint markov
vfoptionsB=struct();
vfoptionsB.exoticpreferences='EpsteinZin';
vfoptionsB.EZutils=0;
vfoptionsB.EZriskaversion='ezgamma';
vfoptionsB.EZeis='ezphi';
simoptionsB=struct();
[VB,PolicyB]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zB,N_j,d_grid,a_grid,z_gridB,pi_zB,ReturnFn_ze_cons,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDistB=StationaryDist_FHorz_Case1(jequaloneDistA,AgeWeightParamNames,PolicyB,n_d,n_a,n_zB,N_j,pi_zB,Params,simoptionsB);

fprintf('Cross test 3: (z,e) vs merged joint-markov, V, should be roughly 1e-13 or smaller: %g \n',max(abs(VA(:)-VB(:))))
fprintf('Cross test 3: (z,e) vs merged joint-markov, Policy, this should be zero: %2.8f \n',max(abs(PolicyA(:)-PolicyB(:))))
fprintf('Cross test 3: (z,e) vs merged joint-markov, StationaryDist, should be roughly zero: %g \n',max(abs(StationaryDistA(:)-StationaryDistB(:))))

clear VA VB PolicyA PolicyB StationaryDistA StationaryDistB

%% Epstein-Zin case 2: Utility-units, positive-valued utility fn (EZutils=1, EZpositiveutility=1)
% Model A: markov z plus iid e
vfoptionsA=struct();
vfoptionsA.exoticpreferences='EpsteinZin';
vfoptionsA.EZutils=1;
vfoptionsA.EZpositiveutility=1;
vfoptionsA.EZriskaversion='ezrisk';
vfoptionsA.n_e=vfoptionsbaseline.n_e;
vfoptionsA.e_grid=vfoptionsbaseline.e_grid;
vfoptionsA.pi_e=vfoptionsbaseline.pi_e;
simoptionsA=struct();
simoptionsA.n_e=simoptionsbaseline.n_e;
simoptionsA.e_grid=simoptionsbaseline.e_grid;
simoptionsA.pi_e=simoptionsbaseline.pi_e;
jequaloneDistA=zeros(n_a,n_z,vfoptionsA.n_e,'gpuArray');
jequaloneDistA(1,ceil(n_z/2),ceil(vfoptionsA.n_e/2))=1; % no assets, midpoint shocks
[VA,PolicyA]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze_posU,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDistA=StationaryDist_FHorz_Case1(jequaloneDistA,AgeWeightParamNames,PolicyA,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

% Model B: the same two shocks as a single joint markov
vfoptionsB=struct();
vfoptionsB.exoticpreferences='EpsteinZin';
vfoptionsB.EZutils=1;
vfoptionsB.EZpositiveutility=1;
vfoptionsB.EZriskaversion='ezrisk';
simoptionsB=struct();
[VB,PolicyB]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zB,N_j,d_grid,a_grid,z_gridB,pi_zB,ReturnFn_ze_posU,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDistB=StationaryDist_FHorz_Case1(jequaloneDistA,AgeWeightParamNames,PolicyB,n_d,n_a,n_zB,N_j,pi_zB,Params,simoptionsB);

fprintf('Cross test 3: (z,e) vs merged joint-markov, V, should be roughly 1e-13 or smaller: %g \n',max(abs(VA(:)-VB(:))))
fprintf('Cross test 3: (z,e) vs merged joint-markov, Policy, this should be zero: %2.8f \n',max(abs(PolicyA(:)-PolicyB(:))))
fprintf('Cross test 3: (z,e) vs merged joint-markov, StationaryDist, should be roughly zero: %g \n',max(abs(StationaryDistA(:)-StationaryDistB(:))))

clear VA VB PolicyA PolicyB StationaryDistA StationaryDistB

%% Epstein-Zin case 3: Utility-units, negative-valued utility fn (EZutils=1, EZpositiveutility=0)
% Model A: markov z plus iid e
vfoptionsA=struct();
vfoptionsA.exoticpreferences='EpsteinZin';
vfoptionsA.EZutils=1;
vfoptionsA.EZpositiveutility=0;
vfoptionsA.EZriskaversion='ezrisk';
vfoptionsA.n_e=vfoptionsbaseline.n_e;
vfoptionsA.e_grid=vfoptionsbaseline.e_grid;
vfoptionsA.pi_e=vfoptionsbaseline.pi_e;
simoptionsA=struct();
simoptionsA.n_e=simoptionsbaseline.n_e;
simoptionsA.e_grid=simoptionsbaseline.e_grid;
simoptionsA.pi_e=simoptionsbaseline.pi_e;
jequaloneDistA=zeros(n_a,n_z,vfoptionsA.n_e,'gpuArray');
jequaloneDistA(1,ceil(n_z/2),ceil(vfoptionsA.n_e/2))=1; % no assets, midpoint shocks
[VA,PolicyA]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze_negU,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDistA=StationaryDist_FHorz_Case1(jequaloneDistA,AgeWeightParamNames,PolicyA,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

% Model B: the same two shocks as a single joint markov
vfoptionsB=struct();
vfoptionsB.exoticpreferences='EpsteinZin';
vfoptionsB.EZutils=1;
vfoptionsB.EZpositiveutility=0;
vfoptionsB.EZriskaversion='ezrisk';
simoptionsB=struct();
[VB,PolicyB]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zB,N_j,d_grid,a_grid,z_gridB,pi_zB,ReturnFn_ze_negU,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDistB=StationaryDist_FHorz_Case1(jequaloneDistA,AgeWeightParamNames,PolicyB,n_d,n_a,n_zB,N_j,pi_zB,Params,simoptionsB);

fprintf('Cross test 3: (z,e) vs merged joint-markov, V, should be roughly 1e-13 or smaller: %g \n',max(abs(VA(:)-VB(:))))
fprintf('Cross test 3: (z,e) vs merged joint-markov, Policy, this should be zero: %2.8f \n',max(abs(PolicyA(:)-PolicyB(:))))
fprintf('Cross test 3: (z,e) vs merged joint-markov, StationaryDist, should be roughly zero: %g \n',max(abs(StationaryDistA(:)-StationaryDistB(:))))

clear VA VB PolicyA PolicyB StationaryDistA StationaryDistB

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
