function output=EstimFHorz_ParametrizeParams(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% EstimFHorzTests part 9: parameters DERIVED by ParametrizeParamsFn feeding the shock process and the initial distribution, through
% EstimateLifeCycleModel_MethodOfMoments. The estimated primitives (half_rho, loginitdecay) appear in NO shock-fn or dist-fn
% signature; the derived rho_z and initdecay do (EstimFHorz_ParametrizeShockDistFn writes them). Before the 2026-10-12 fix the
% calibrateshocks detection name-matched only EstimParamNames, so the shocks were SILENTLY FROZEN at the initial parameters here;
% the fix matches against EstimParamNames plus the ParametrizeParamsFn outputs (estimoptions.ParametrizeParamsFn_OutputNames,
% filled by ParseParametrizeParamsFn when left empty). The grid interpolation layer is on.
%   (a) ParseParametrizeParamsFn directly: parsefailed=0, the discovered inputs are exactly {half_rho, loginitdecay}, the outputs
%       contain the derived and the input names; a function reading a non-existent parameter gives parsefailed=1
%   (b) vfoptions.ExogShockFn(rho_z,sigma_z) with rho_z=2*half_rho derived and half_rho estimated (a STALE rho_z sits in Params to
%       prove it is re-derived): at the truth (skipestimation=1) the objective is zero and the Jacobian equals finite differences
%       the test takes itself with ParametrizeParamsFn applied and the grids rebuilt (exact); the silent-freeze regression check
%       (J numerically zero would be the old bug); user-set estimoptions.calibrateshocks=0 is honored (J zero, the shocks stay at
%       the stale rho_z) and fires a warning; user-given estimoptions.ParametrizeParamsFn_OutputNames omitting rho_z is used
%       verbatim (shocks frozen, J zero)
%   (c) jequaloneDist as a function taking the derived initdecay=exp(loginitdecay), loginitdecay estimated (stale initdecay in
%       Params): at the truth the objective is zero (pins that ParametrizeParamsFn runs before jequaloneDist is evaluated) and J
%       equals finite differences the test takes itself (exact); started away it is recovered (close)
%   (d) half_rho started away from the truth is recovered (close)
% The shock function is the CalibFHorzTests one (on the path). Each driver run is wrapped so an error in the plumbing is reported
% as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
Jr=Params.Jr;
vfoptions=struct(); simoptions=struct();
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
simoptions_full=simoptions; simoptions_full.whichstats=ones(1,7);
estimoptions=struct(); estimoptions.verbose=0;
ParametrizeParamsFn=@(P) EstimFHorz_ParametrizeShockDistFn(P);

%% (a) ParseParametrizeParamsFn directly
Params_p=Params; Params_p.half_rho=0.45; Params_p.loginitdecay=log(0.7);
[inN,outN,pfail]=ParseParametrizeParamsFn(ParametrizeParamsFn,Params_p);
d1=double(~isequal(sort(inN),sort({'half_rho','loginitdecay'})));
d2=double(~all(ismember({'rho_z','initdecay','half_rho','loginitdecay'},outN)));
fprintf('ParametrizeParams (a) ParseParametrizeParamsFn: parsefailed, this should be zero: %.3e; inputs are exactly half_rho and loginitdecay, this should be zero: %.3e; outputs contain the derived and the input names, this should be zero: %.3e \n',pfail,d1,d2)
[~,~,pfail2]=ParseParametrizeParamsFn(@(P) setfield(P,'x',P.thisfielddoesnotexist),Params_p);
fprintf('ParametrizeParams (a) a fn reading a non-existent parameter: 1-parsefailed, this should be zero: %.3e \n',1-pfail2)

%% (b) ExogShockFn fed by the derived rho_z, half_rho estimated
% Truth: rho_z=2*0.45=0.9, sigma_z=0.03 -- the bank's own z grid is the Farmer-Toda discretisation at exactly these values, so the
% targets can be built from a plain solve on the bank grids
Params_b=Params_p; Params_b.sigma_z=0.03;
Params_b.rho_z=0.123; % STALE on purpose: the command must re-derive rho_z=0.9 from half_rho via ParametrizeParamsFn before the shocks are set up
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_b,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params_b,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params_b,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params_b,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
TMb=struct();
t=gather(LCPFull.earnings.StdDeviation); t(Jr:end)=NaN; TMb.AgeConditionalStats.earnings.StdDeviation=t; % (earnings are zero in retirement)
TMb.AllStats.earnings.Gini=gather(AllFull.earnings.Gini);
nT=(Jr-1)+1; W=eye(nT); Omega=diag(0.01*(1+(1:nT)'/nT));
vfoptions_b=vfoptions; vfoptions_b.ExogShockFn=@(rho_z,sigma_z) CalibFHorz_ExogShockFn(rho_z,sigma_z);
simoptions_b=simoptions; simoptions_b.ExogShockFn=vfoptions_b.ExogShockFn;
estimoptions_s=estimoptions; estimoptions_s.skipestimation=1;
estimoptions_s.eedefault=3; % the checks below are at the epsilon10000 step (the default is 'auto')
ran=0; objv=NaN; maxd=NaN; freezefail=NaN;
try
    [~,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments({'half_rho'},TMb,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_b,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,ParametrizeParamsFn,FnsToEvaluate,estimoptions_s,vfoptions_b,simoptions_b);
    ran=1; objv=estsummary.objectivefnval;
    dc=estsummary.doublechecks.epsilon10000;
    % finite differences by hand: ParametrizeParamsFn applied, shocks rebuilt (forced: the test knows the shocks depend on half_rho)
    mu=EstimFHorz_ModelMoments(dc.estimparamsvecup(:),{'half_rho'},TMb,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_b,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,FnsToEvaluate,vfoptions_b,simoptions_b,0,[],ParametrizeParamsFn,1);
    md=EstimFHorz_ModelMoments(dc.estimparamsvecdown(:),{'half_rho'},TMb,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_b,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,FnsToEvaluate,vfoptions_b,simoptions_b,0,[],ParametrizeParamsFn,1);
    Jc=(mu-md)/(dc.estimparamsvecup(1)-dc.estimparamsvecdown(1));
    maxd=CoreMaxAbs(dc.J_centered(:)-Jc(:));
    freezefail=double(CoreMaxAbs(dc.J_centered)<1e-8); % the old bug: shocks frozen, so the moments do not respond to half_rho at all
catch err
    fprintf('ParametrizeParams (b): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ParametrizeParams (b) derived rho_z, half_rho, skipestimation at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; J_centered vs finite differences with ParametrizeParamsFn and the grids rebuilt, this should be zero: %.3e; silent shock freeze (J numerically zero), this should be zero: %.3e \n',1-ran,objv,maxd,freezefail)
% user-set calibrateshocks=0 is honored (shocks stay at the stale rho_z, J is exactly zero) and fires a warning. The warning is
% read from the captured command-window output (evalc), NOT from lastwarn: the frozen J=0 makes the standard-error sandwich invert
% a singular matrix, whose own warning fires later and overwrites lastwarn (that cost a bank run, 2026-10-12).
estimoptions_o=estimoptions_s; estimoptions_o.calibrateshocks=0; %#ok<STRNU> (used inside the evalc string below, which static analysis cannot see)
ran=0; dJ0=NaN; warned=NaN;
try
    lastwarn('');
    cap=evalc('[~,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments({''half_rho''},TMb,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_b,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,ParametrizeParamsFn,FnsToEvaluate,estimoptions_o,vfoptions_b,simoptions_b);');
    ran=1; dJ0=CoreMaxAbs(estsummary.doublechecks.epsilon10000.J_centered);
    warned=double(contains(cap,'calibrateshocks')||contains(lastwarn,'calibrateshocks')); % the contradicted-override warning, wherever it surfaced
catch err
    fprintf('ParametrizeParams (b) calibrateshocks=0: the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ParametrizeParams (b) user-set calibrateshocks=0: did not run, this should be zero: %.3e; J with the shocks frozen as instructed, this should be zero: %.3e; missing contradicted-override warning, this should be zero: %.3e \n',1-ran,dJ0,1-warned)
% user-given ParametrizeParamsFn_OutputNames omitting rho_z is used verbatim (no re-parse): the shocks freeze
estimoptions_u=estimoptions_s; estimoptions_u.ParametrizeParamsFn_OutputNames={'initdecay'}; estimoptions_u.ParametrizeParamsFn_InputNames={'half_rho','loginitdecay'};
ran=0; dJu=NaN;
try
    [~,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments({'half_rho'},TMb,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_b,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,ParametrizeParamsFn,FnsToEvaluate,estimoptions_u,vfoptions_b,simoptions_b);
    ran=1; dJu=CoreMaxAbs(estsummary.doublechecks.epsilon10000.J_centered);
catch err
    fprintf('ParametrizeParams (b) user _OutputNames: the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ParametrizeParams (b) user-given _OutputNames omitting rho_z (used verbatim, so the shocks freeze): did not run, this should be zero: %.3e; J, this should be zero: %.3e \n',1-ran,dJu)

%% (d) recovery of half_rho from a wrong start (0.40, truth 0.45)
Params_start=Params_b; Params_start.half_rho=0.40;
estimoptions_r=estimoptions; estimoptions_r.fminalgo=8; estimoptions_r.constrain0to1={'half_rho'};
ran=0; perr=NaN; tt=NaN; objv=NaN;
try
    t0=tic;
    [EstimParams,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments({'half_rho'},TMb,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,ParametrizeParamsFn,FnsToEvaluate,estimoptions_r,vfoptions_b,simoptions_b);
    tt=toc(t0); ran=1; perr=abs(EstimParams.half_rho-0.45); objv=estsummary.objectivefnval;
catch err
    fprintf('ParametrizeParams (d) from 0.40: the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ParametrizeParams (d) half_rho started at 0.40 (truth 0.45) (%.1f s), objective %.3e: did not run, this should be zero: %.3e; half_rho error, this should be close to zero: %.3e \n',tt,objv,1-ran,perr)

%% (c) jequaloneDist as a function taking the derived initdecay, loginitdecay estimated
jequaloneDistFn=@(a_grid,z_grid,n_a,n_z,initdecay) EstimFHorz_jequaloneDistParamFn(a_grid,z_grid,n_a,n_z,initdecay);
Params_c=Params_p; Params_c.loginitdecay=log(0.9); % truth: initdecay=0.9 (geometric over the first 60 asset points, see EstimFHorz_jequaloneDistParamFn on why identification needs the wide support)
Params_c.initdecay=0.123; % STALE on purpose: the objective must re-derive initdecay=0.9 before evaluating jequaloneDist
jequaloneDist_truth=gpuArray(EstimFHorz_jequaloneDistParamFn(a_grid,z_grid,n_a,n_z,0.9));
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_c,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist_truth,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params_c,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params_c,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
TMc=struct();
t=gather(LCPFull.assets.Mean); t(t<1e-6)=NaN; TMc.AgeConditionalStats.assets.Mean=t; % age 1 stays in: its mean assets respond to initdecay directly
TMc.AllStats.assets.Gini=gather(AllFull.assets.Gini);
nTc=sum(~isnan(t))+1; Wc=eye(nTc); Omegac=diag(0.01*(1+(1:nTc)'/nTc));
simoptions_c=simoptions; simoptions_c.jequaloneDist_usergrids=1;
simoptions_c.a_grid=a_grid; simoptions_c.z_grid=z_grid; % (a functional jequaloneDist needs the grids in simoptions)
ran=0; objv=NaN; maxd=NaN; flatfail=NaN;
try
    [~,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments({'loginitdecay'},TMc,Wc,Omegac,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_c,DiscountFactorParamNames,jequaloneDistFn,AgeWeightParamNames,ParametrizeParamsFn,FnsToEvaluate,estimoptions_s,vfoptions,simoptions_c);
    ran=1; objv=estsummary.objectivefnval;
    dc=estsummary.doublechecks.epsilon10000;
    mu=EstimFHorz_ModelMoments(dc.estimparamsvecup(:),{'loginitdecay'},TMc,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_c,DiscountFactorParamNames,jequaloneDistFn,AgeWeightParamNames,FnsToEvaluate,vfoptions,simoptions_c,0,[],ParametrizeParamsFn);
    md=EstimFHorz_ModelMoments(dc.estimparamsvecdown(:),{'loginitdecay'},TMc,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_c,DiscountFactorParamNames,jequaloneDistFn,AgeWeightParamNames,FnsToEvaluate,vfoptions,simoptions_c,0,[],ParametrizeParamsFn);
    Jc=(mu-md)/(dc.estimparamsvecup(1)-dc.estimparamsvecdown(1));
    maxd=CoreMaxAbs(dc.J_centered(:)-Jc(:));
    flatfail=double(CoreMaxAbs(dc.J_centered)<1e-6); % guards against a frozen initdecay AND against this test degenerating into an unidentified (flat) one, as its first version did
catch err
    fprintf('ParametrizeParams (c): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ParametrizeParams (c) jequaloneDist fn taking the derived initdecay, skipestimation at the truth: did not run, this should be zero: %.3e; objective (pins ParametrizeParamsFn before jequaloneDist), this should be zero: %.3e; J_centered vs finite differences, this should be zero: %.3e; flat objective (frozen or unidentified), this should be zero: %.3e \n',1-ran,objv,maxd,flatfail)
Params_start=Params_c; Params_start.loginitdecay=log(0.8);
estimoptions_r2=estimoptions; estimoptions_r2.fminalgo=8;
ran=0; perr=NaN; tt=NaN; objv=NaN;
try
    t0=tic;
    [EstimParams,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments({'loginitdecay'},TMc,Wc,Omegac,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDistFn,AgeWeightParamNames,ParametrizeParamsFn,FnsToEvaluate,estimoptions_r2,vfoptions,simoptions_c);
    tt=toc(t0); ran=1; perr=abs(EstimParams.loginitdecay-log(0.9)); objv=estsummary.objectivefnval;
catch err
    fprintf('ParametrizeParams (c) from log(0.8): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ParametrizeParams (c) loginitdecay started at log(0.8) (truth log(0.9)) (%.1f s), objective %.3e: did not run, this should be zero: %.3e; loginitdecay error, this should be close to zero: %.3e \n',tt,objv,1-ran,perr)

output=struct();
output.ran=ran;

end
