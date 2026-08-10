function output=CoreInfHorzVFIAlgo_CrossTest_dummyd(n_a,n_a_big,n_z,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)
% CROSS-TEST: a decision variable d that DOES NOTHING vs no d at all.
%
% Model A: has a d variable with n_d_A grid points that are ALL the same value (dval),
%          so the choice of d is irrelevant -- economically it is a no-op.
% Model B: has NO d, but a parameter dparam=dval takes d's place in the (same) return fn.
%
% The two models are literally the same economic problem, so V must be identical (and the
% aprime policy identical). Model A forces the with-d 'Refine' code path; Model B uses the
% nod code path. So this isolates the Refine path: any difference -- especially one that only
% shows up WITH grid interpolation -- points at the Refine/GI code rather than economics.
%
% Uses the z (markov), no-e model.

DF=DiscountFactorParamNames;

dval=0.5; % interior value (need d<1 for the leisure term); all d grid points equal this
Params.dparam=dval;

% Model A: dummy d (5 grid points, all equal to dval)
n_d_A=5;
d_grid_A=dval*ones(n_d_A,1);
ReturnFn_A=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);

% Model B: no d; dparam plays the role of d in the SAME return function
ReturnFn_B=@(aprime,a,z,r,w,sigma,eta,varphi,dparam) ReturnFn_d_z_noe_nosemiz(dparam,aprime,a,z,r,w,sigma,eta,varphi);

fprintf('\n================ CROSS-TEST: dummy-d (n_d=%d, all d=%g) vs no-d (dparam=%g) ================\n',n_d_A,dval,dval);

%% Without GI
vfo=struct(); vfo.verbose_advice=0; % silence the postGI maxaprimediff advice in this cross-test
[VA,PolicyA]=ValueFnIter_InfHorz(n_d_A,n_a,n_z,d_grid_A,a_grid,z_grid,pi_z,ReturnFn_A,Params,DF,[],vfo);
[VB,PolicyB]=ValueFnIter_InfHorz(0,    n_a,n_z,[],      a_grid,z_grid,pi_z,ReturnFn_B,Params,DF,[],vfo);
% PolicyA is [d;aprime] (2 rows); PolicyB is [aprime] (1 row). Compare V, and the aprime index.
aprimeA=PolicyA(2,:,:);
fprintf('noGI: dummy-d vs no-d, V      should be ~0: %2.8f \n',max(abs(VA(:)-VB(:))));
fprintf('noGI: dummy-d vs no-d, aprime should be  0: %2.8f \n',max(abs(aprimeA(:)-PolicyB(:))));

%% With GI (postGI, the default)
vfo=struct(); vfo.gridinterplayer=1; vfo.ngridinterp=5; vfo.verbose_advice=0;
[VAg,PolicyAg]=ValueFnIter_InfHorz(n_d_A,n_a,n_z,d_grid_A,a_grid,z_grid,pi_z,ReturnFn_A,Params,DF,[],vfo);
[VBg,PolicyBg]=ValueFnIter_InfHorz(0,    n_a,n_z,[],      a_grid,z_grid,pi_z,ReturnFn_B,Params,DF,[],vfo);
fprintf('GI  (postGI): dummy-d vs no-d, V should be ~0: %2.8f \n',max(abs(VAg(:)-VBg(:))));

%% With GI (preGI)
vfo=struct(); vfo.gridinterplayer=1; vfo.ngridinterp=5; vfo.preGI=1; vfo.verbose_advice=0;
[VApre,~]=ValueFnIter_InfHorz(n_d_A,n_a,n_z,d_grid_A,a_grid,z_grid,pi_z,ReturnFn_A,Params,DF,[],vfo);
[VBpre,~]=ValueFnIter_InfHorz(0,    n_a,n_z,[],      a_grid,z_grid,pi_z,ReturnFn_B,Params,DF,[],vfo);
fprintf('GI  (preGI):  dummy-d vs no-d, V should be ~0: %2.8f \n',max(abs(VApre(:)-VBpre(:))));

output=struct();

end
