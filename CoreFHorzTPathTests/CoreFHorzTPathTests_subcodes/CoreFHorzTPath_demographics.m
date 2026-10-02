function output=CoreFHorzTPath_demographics(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,AgeWeightParamNames,transpathoptionsbaseline)
% Checks that changing population demographics are handled correctly on an FHorz transition path:
% age weights (mewj) and conditional survival (sj) that change over time, and the distribution of
% newborns (jequalOneDist) as a path. On the life-cycle model of CoreFHorzTPath_nod_z_noe_nosemiz.
%
% How the FHorz path treats them (TransitionPath_Case1_FHorz):
%   The agent distribution is normalised to one within each age, and the age weights only enter when
%     aggregating, so an aggregate in period t is sum_j mewj(j,t)*E[f|age j, period t]. When mewj is on the
%     ParamPath (N_j-by-T) it is read period by period; the user supplies it, nothing derives it from sj.
%   sj only enters the value fn, as one of the DiscountFactorParamNames, read period by period.
%   jequalOneDist can be a T-path; the cohort that is age 1 in period t is jequalOneDist(:,t) (t>=2,
%     period 1 is AgentDist_initial).
%
% THE METHOD is that of doPart(6): exact oracles, the path evaluated without moving it (shooting with
% maxiter=1 returns the conditions at the path it was handed), and an aggregate exposed as the path's
% general eqm condition by a GE eqn that just returns it.
%   (1) mewj on the ParamPath vs mewj held at period 1 with the ratio mewj(j,t)/mewj(j,1) folded into the
%       FnsToEvaluate as an age-dependent ParamPath parameter: the same sum, so the same aggregates
%   (2) mewj on the ParamPath but constant vs not on the ParamPath (ageweightstrivial=0 vs =1)
%   (3) (1) and (2) with fastOLG=1 and 0
%   (4) sj on the ParamPath in DiscountFactorParamNames={'beta','sj'} vs one discount factor betasj=beta*sj
%   (5) sj on the ParamPath but constant vs not on the ParamPath
%   (6) jequalOneDist as a T-path: the age-1 aggregate of assets in period t is mewj(1)*a_grid(newborn
%       asset point in period t), which can be written down by hand; this pins the timing
%   (7) jequalOneDist constant but given as a T-path vs given once
%   (8) AgentDistOnTransPath + EvalFnOnTransPath_AggVars with the mewj path vs the aggregates inside the path
%   (9) a general eqm path through the whole demographic transition (sj and mewj both changing), with
%       bequests (1-sj)*aprime received next period: Anderson vs shooting, and Anderson's answer
%       re-evaluated with the bequest path shifted by hand

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];
DiscountFactorParamNames={'beta'};

% zero assets, mid point for the shock
jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;

%% The demographic transition: survival rises and births slow down, over the first half of the path
sj_initial=[linspace(0.99,0.9,N_j-1),0]';
sj_final=[linspace(0.995,0.95,N_j-1),0]';
n_initial=0.01;
n_final=0;
ramp=min((0:T-1)/(T/2),1); % 0 to 1 over the first half of the path
sjPath=sj_initial+(sj_final-sj_initial).*ramp; % N_j-by-T
nPath=n_initial+(n_final-n_initial)*ramp; % 1-by-T, growth of the newborn cohort
% Mass by age. Period 1 is the stationary population of the initial demographics. After that each age
% is last period's previous age times its survival, and the newborn cohort grows at n.
massPath=zeros(N_j,T);
massPath(1,1)=1;
for jj=2:N_j
    massPath(jj,1)=massPath(jj-1,1)*sj_initial(jj-1)/(1+n_initial);
end
for tt=2:T
    massPath(1,tt)=massPath(1,tt-1)*(1+nPath(tt));
    massPath(2:N_j,tt)=massPath(1:N_j-1,tt-1).*sjPath(1:N_j-1,tt-1);
end
mewjPath=massPath./sum(massPath,1); % N_j-by-T
fprintf('Demographics: the age weights move by up to %.3e over the path (must not be zero, or (1) is empty) \n',max(max(abs(mewjPath-mewjPath(:,1)))))

% Period-1 demographics in Params, so that the initial agent distribution agrees with period 1 of the path
Params.mewj=mewjPath(:,1)';
Params.sj=sj_initial';

[V_final,Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_final,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

% Evaluate-only: one shooting iteration returns the conditions at the path it was given
transpathoptionsE=transpathoptionsbaseline;
transpathoptionsE.maxiter=1;
transpathoptionsE.GEnewprice=3;
transpathoptionsE.GEnewprice3.howtoupdate={'g','r',0,0.1};
transpathoptionsE.verbose=0;

GeneralEqmEqnsAE.g=@(assets,earnings) assets+10*earnings; % exposes both aggregates at once (earnings carries the age-dependent kappa_j)

%% (1)-(3) Age weights on the ParamPath
fastOLGopts=[1,0];
for ff=1:length(fastOLGopts)
    transpathoptionsE.fastOLG=fastOLGopts(ff);

    % (1) mewj on the ParamPath
    ParamPath1=ParamPath;
    ParamPath1.mewj=mewjPath;
    [~,GEcondnPath1]=TransitionPath_Case1_FHorz(PricePath, ParamPath1, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsAE, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    % Oracle: mewj held at period 1, with mewj(j,t)/mewj(j,1) in the FnsToEvaluate
    ParamPath1O=ParamPath;
    ParamPath1O.mewjratio=mewjPath./mewjPath(:,1); % N_j-by-T, an age-dependent parameter on the path
    Params1O=Params;
    Params1O.mewjratio=ones(1,N_j);
    FnsToEvaluate1O.assets=@(aprime,a,z,mewjratio) a*mewjratio;
    FnsToEvaluate1O.earnings=@(aprime,a,z,w,kappa_j,mewjratio) w*kappa_j*z*mewjratio;
    [~,GEcondnPath1O]=TransitionPath_Case1_FHorz(PricePath, ParamPath1O, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate1O, GeneralEqmEqnsAE, Params1O, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('Demographics (fastOLG=%i): (1) mewj on the ParamPath vs period-1 mewj with the ratio in the FnsToEvaluate, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(GEcondnPath1.g-GEcondnPath1O.g)))
    if ff==1
        GEcondnPath1keep=GEcondnPath1; % for (8)
    end

    % (2) mewj on the ParamPath but constant, vs not on the ParamPath at all
    ParamPath2=ParamPath;
    ParamPath2.mewj=repmat(mewjPath(:,1),1,T);
    [~,GEcondnPath2]=TransitionPath_Case1_FHorz(PricePath, ParamPath2, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsAE, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    [~,GEcondnPath2O]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsAE, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('Demographics (fastOLG=%i): (2) constant mewj on the ParamPath vs not on it, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(GEcondnPath2.g-GEcondnPath2O.g)))
end
transpathoptionsE.fastOLG=transpathoptionsbaseline.fastOLG;

%% (4)-(5) Survival on the ParamPath, in the discount factor
% Same value fn either way: the discount factor is the product of the DiscountFactorParamNames.
transpathoptionsV=transpathoptionsbaseline;
ParamPath4=ParamPath;
ParamPath4.sj=sjPath;
[VPath4,PolicyPath4]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath4, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, {'beta','sj'}, ReturnFn, transpathoptionsV, vfoptions);
ParamPath4O=ParamPath;
ParamPath4O.betasj=Params.beta*sjPath;
Params4O=Params;
Params4O.betasj=Params.beta*Params.sj;
[VPath4O,PolicyPath4O]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath4O, T, V_final, Policy_final, Params4O, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, {'betasj'}, ReturnFn, transpathoptionsV, vfoptions);
fprintf('Demographics: (4) sj on the ParamPath in the discount factor vs betasj=beta*sj, Policy, this should be zero: %.3e \n',max(abs(PolicyPath4(:)-PolicyPath4O(:))))
fprintf('Demographics: (4) sj on the ParamPath in the discount factor vs betasj=beta*sj, V, this should be zero: %.3e \n',max(abs(VPath4(:)-VPath4O(:))))
% and the sj path must matter, or the comparison above is empty
[VPath4C,~]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, {'beta','sj'}, ReturnFn, transpathoptionsV, vfoptions);
fprintf('Demographics: for scale, V with the sj path vs with sj held at period 1 differs by %.3e (must not be zero, or (4) is empty) \n',max(abs(VPath4(:)-VPath4C(:))))
% (5) constant sj on the ParamPath vs not on it
ParamPath5=ParamPath;
ParamPath5.sj=repmat(sjPath(:,1),1,T);
[VPath5,PolicyPath5]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath5, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, {'beta','sj'}, ReturnFn, transpathoptionsV, vfoptions);
[VPath5O,PolicyPath5O]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, {'beta','sj'}, ReturnFn, transpathoptionsV, vfoptions);
fprintf('Demographics: (5) constant sj on the ParamPath vs not on it, V and Policy, this should be zero: V %.3e, Policy %.3e \n',max(abs(VPath5(:)-VPath5O(:))),max(abs(PolicyPath5(:)-PolicyPath5O(:))))

%% (6)-(7) The distribution of newborns as a T-path
% Newborns enter at asset grid point ia_t, which climbs over the path. The age-1 aggregate of assets in
% period t>=2 is then mewj(1)*a_grid(ia_t), exactly. Period 1 is AgentDist_initial, so it is not compared.
ia_T=1+floor(20*(0:T-1)/T); % grid points 1 to 20
jequaloneDist_T=zeros([n_a,n_z,T],'gpuArray');
for tt=1:T
    jequaloneDist_T(ia_T(tt),ceil(n_z/2),tt)=1;
end
FnsToEvaluate6.assets1=@(aprime,a,z,agej) a*(agej==1);
GeneralEqmEqns6.g=@(assets1) assets1;
[~,GEcondnPath6]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist_T, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate6, GeneralEqmEqns6, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
Assets1byhand=Params.mewj(1)*gather(a_grid(ia_T(2:T-1)))';
fprintf('Demographics: (6) jequaloneDist as a T-path, age-1 assets in periods 2 to T-1 vs by hand, this should be zero: %.3e \n',max(abs(gather(GEcondnPath6.g(2:T-1))-Assets1byhand)))
fprintf('Demographics: for scale, age-1 assets move by %.3e over the path (must not be zero, or (6) is empty) \n',max(Assets1byhand)-min(Assets1byhand))
% (7) constant newborns given as a T-path vs given once
[~,GEcondnPath7]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, repmat(jequaloneDist,1,1,T), n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsAE, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
[~,GEcondnPath7O]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsAE, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
fprintf('Demographics: (7) constant jequaloneDist as a T-path vs given once, this should be zero: %.3e \n',max(abs(GEcondnPath7.g-GEcondnPath7O.g)))

%% (8) The standalone commands, with the mewj path
transpathoptions8=transpathoptionsbaseline;
[~,PolicyPath8]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath1, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions8, vfoptions);
AgentDistPath8=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPath1, PolicyPath8, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptions8, simoptions);
AggVarsPath8=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath8, PolicyPath8, PricePath, ParamPath1, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptions8, simoptions);
fprintf('Demographics: (8) AgentDistOnTransPath+EvalFnOnTransPath_AggVars with the mewj path vs inside the path, this should be zero: %.3e \n',max(abs(gather(AggVarsPath8.assets.Mean(1:T-1)+10*AggVarsPath8.earnings.Mean(1:T-1))-gather(GEcondnPath1keep.g))))

%% (12) The agent dist laid out the slow way (simoptions.fastOLG=0), inside the path
% It is then (a,z,j) rather than (a,j,z), and Step4tt_AggVars puts it into the (a,j,z) order the
% FnsToEvaluate are evaluated in before aggregating. Before that, nothing did, and nothing tested it.
transpathoptions12=transpathoptionsE;
transpathoptions12.fastOLG=0; % simoptions.fastOLG=0 requires transpathoptions.fastOLG=0
simoptions12=simoptions;
simoptions12.fastOLG=0;
[~,GEcondnPath12]=TransitionPath_Case1_FHorz(PricePath, ParamPath1, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsAE, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions12, simoptions12, vfoptions);
fprintf('Demographics: (12) simoptions.fastOLG=0 vs =1, aggregates inside the path with the mewj path, this should be zero: %.3e \n',max(abs(GEcondnPath12.g-GEcondnPath1keep.g)))

%% (9) A general eqm path through the demographic transition
% r_t = 0.05*(1+K_t)^(-0.5), with K_t assets plus last period's bequests, and bequests (1-sj)*aprime.
% sj is on the ParamPath twice over: in the discount factor and in the bequests.
FnsToEvaluate9=FnsToEvaluate;
FnsToEvaluate9.beq=@(aprime,a,z,sj) (1-sj)*aprime;
GeneralEqmEqns9.rmarket=@(r,assets,beq_tminus1) r-0.05*(1+assets+beq_tminus1)^(-0.5);
ParamPath9=ParamPath;
ParamPath9.mewj=mewjPath;
ParamPath9.sj=sjPath;
initialvalues9.beq=0.02;

transpathoptions9A=transpathoptionsbaseline;
transpathoptions9A.GEnewprice=2;
transpathoptions9A.GEnewprice2.howtoupdate={'rmarket','r',0,0.5};
transpathoptions9A.maxiter=200;
transpathoptions9A.initialvalues=initialvalues9;
transpathoptions9A.verbose=0;
[PricePath9A,GEcondnPath9A]=TransitionPath_Case1_FHorz(PricePath, ParamPath9, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate9, GeneralEqmEqns9, Params, {'beta','sj'}, AgeWeightParamNames, transpathoptions9A, simoptions, vfoptions);
transpathoptions9B=rmfield(transpathoptions9A,'GEnewprice2');
transpathoptions9B.GEnewprice=3;
transpathoptions9B.GEnewprice3.howtoupdate={'rmarket','r',0,0.5};
transpathoptions9B.maxiter=500;
[PricePath9B,GEcondnPath9B]=TransitionPath_Case1_FHorz(PricePath, ParamPath9, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate9, GeneralEqmEqns9, Params, {'beta','sj'}, AgeWeightParamNames, transpathoptions9B, simoptions, vfoptions);
fprintf('Demographics GE path: how far Anderson moved r from the starting guess: %.3e \n',max(abs(PricePath9A.r-PricePath.r)))
fprintf('Demographics GE path, Anderson, GE condn distance: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',max(abs(GEcondnPath9A.rmarket)),1e-4)
fprintf('Demographics GE path, shooting, GE condn distance: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',max(abs(GEcondnPath9B.rmarket)),1e-4)
fprintf('Demographics GE path, Anderson vs shooting, same eqm path, should be small: %.3e \n',max(abs(PricePath9A.r-PricePath9B.r)))
% Anderson's answer, re-evaluated with the bequest path shifted by hand
transpathoptions9E=transpathoptionsE;
GeneralEqmEqns9beq.g=@(beq) beq;
[~,GEcondnPath9beq]=TransitionPath_Case1_FHorz(PricePath9A, ParamPath9, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate9, GeneralEqmEqns9beq, Params, {'beta','sj'}, AgeWeightParamNames, transpathoptions9E, simoptions, vfoptions);
beqPath=gather(GEcondnPath9beq.g); % periods 1,...,T-1
ParamPath9O=ParamPath9;
ParamPath9O.beqlag=[initialvalues9.beq, beqPath(1:T-2), beqPath(T-1)]; % the last entry is never read
Params9O=Params;
Params9O.beqlag=initialvalues9.beq;
GeneralEqmEqns9O.g=@(r,assets,beqlag) r-0.05*(1+assets+beqlag)^(-0.5);
[~,GEcondnPath9O]=TransitionPath_Case1_FHorz(PricePath9A, ParamPath9O, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate9, GeneralEqmEqns9O, Params9O, {'beta','sj'}, AgeWeightParamNames, transpathoptions9E, simoptions, vfoptions);
fprintf('Demographics GE path: conditions Anderson returned vs its answer re-evaluated with the bequests shifted by hand, this should be zero: %.3e \n',max(abs(GEcondnPath9A.rmarket-GEcondnPath9O.g)))

output=struct();

end
