function output=CoreFHorzTPathExpAssetz_timevaryinggrids(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline)
% z and e grids that vary over the transition path, with an experience asset whose aprimeFn depends on z.
%
% THE METHOD. Two models that are the same model:
%   (A) z_grid as [n_z,T] and e_grid as [n_e,T], drifting up by a factor 1+0.1*t/T
%   (B) the plain z_grid and e_grid, with the same drift entering as parameters zscale and escale on ParamPath:
%       the ReturnFn, the aprimeFn and the FnsToEvaluate all see z*zscale and e*escale
% In both, pi_z and pi_e move over the path (a third of the way towards uniform by period T). So everything
% must agree exactly: (A) is the shock setup's time-varying grids, (B) does not use them at all. In
% particular the aprimeFn in (A) must be evaluated on each period's z grid, and in (B) on each period's
% zscale, in the value fn, the agent dist and the GE path alike.
%   (X1) value fn path, (A) vs (B), fastOLG=0 and fastOLG=1
%   (X2) agent dist path, (A) vs (B), fastOLG=0 and fastOLG=1
%   (X3) EvalFnOnTransPath_AggVars, (A) vs (B)
%   (X4) one iteration of the GE path (a general eqm condition that depends on earnings), (A) vs (B),
%        transpathoptions.fastOLG=0 and 1 (the agent dist stays at fastOLG=1, as the GE path requires here)
%   and, for scale, how far the drift moves V and how far one GE iteration moves r (must not be zero)

n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;
driftfactor=1+0.1*(1:1:T)/T; % [1,T]
z_gridT=z_grid.*driftfactor; % [n_z,T]
e_gridT=e_grid.*driftfactor; % [n_e,T]
pimixT=0.3*(1:1:T)/T; % [1,T]
pi_zT=(1-reshape(pimixT,[1,1,T])).*pi_z+reshape(pimixT,[1,1,T]).*ones(n_z,n_z)/n_z; % [n_z,n_z,T]
pi_eT=(1-pimixT).*pi_e+pimixT.*ones(n_e,1)/n_e; % [n_e,T]

%% (A) the drifting grids
vfoptionsA=struct();
vfoptionsA.experienceassetz=vfoptionsbaseline.experienceassetz;
vfoptionsA.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptionsA.n_e=n_e; vfoptionsA.e_grid=e_gridT; vfoptionsA.pi_e=pi_eT;
simoptionsA=struct();
simoptionsA.experienceassetz=simoptionsbaseline.experienceassetz;
simoptionsA.aprimeFn=simoptionsbaseline.aprimeFn;
simoptionsA.d_grid=d_grid;
simoptionsA.a_grid=a_grid;
simoptionsA.z_grid=z_gridT; % aprimeFn depends on z, so sim codes need z_grid
simoptionsA.n_e=n_e; simoptionsA.e_grid=e_gridT; simoptionsA.pi_e=pi_eT;
ReturnFnA=@(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_TPathExpAssetz_nod1_z_e_nosemiz(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluateA=struct();
FnsToEvaluateA.assets=@(d2,a1prime,a1,a2,z,e) a1+a2;
FnsToEvaluateA.earnings=@(d2,a1prime,a1,a2,z,e,w,kappa_j) w*kappa_j*d2*a2*z*e;

%% (B) the plain grids, with the drift as parameters on the path
vfoptionsB=struct();
vfoptionsB.experienceassetz=vfoptionsbaseline.experienceassetz;
vfoptionsB.aprimeFn=@(d2,a2,z,phi1,phi2,zscale) aprimeFn_CoreTestTPathExpAssetz(d2,a2,z*zscale,phi1,phi2);
vfoptionsB.n_e=n_e; vfoptionsB.e_grid=e_grid; vfoptionsB.pi_e=pi_eT;
simoptionsB=struct();
simoptionsB.experienceassetz=simoptionsbaseline.experienceassetz;
simoptionsB.aprimeFn=vfoptionsB.aprimeFn;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;
simoptionsB.z_grid=z_grid;
simoptionsB.n_e=n_e; simoptionsB.e_grid=e_grid; simoptionsB.pi_e=pi_eT;
ReturnFnB=@(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension,zscale,escale) ReturnFn_TPathExpAssetz_nod1_z_e_nosemiz(d2,a1prime,a1,a2,z*zscale,e*escale,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluateB=struct();
FnsToEvaluateB.assets=@(d2,a1prime,a1,a2,z,e) a1+a2;
FnsToEvaluateB.earnings=@(d2,a1prime,a1,a2,z,e,w,kappa_j,zscale,escale) w*kappa_j*d2*a2*z*zscale*e*escale;
ParamPathB=ParamPath;
ParamPathB.zscale=driftfactor;
ParamPathB.escale=driftfactor;
ParamsB=Params;
ParamsB.zscale=1;
ParamsB.escale=1;

%% The final period, and the initial dist (the same for both, from the plain model)
V_final=zeros([n_a,n_z,n_e,N_j],'gpuArray');
Policy_final=ones([2,n_a,n_z,n_e,N_j],'gpuArray');
jequaloneDist=zeros([n_a,n_z,n_e],'gpuArray');
jequaloneDist(1,1,ceil(n_z/2),ceil(n_e/2))=1;
simoptions0=simoptionsA;
simoptions0.z_grid=z_grid; simoptions0.e_grid=e_grid; simoptions0.pi_e=pi_e;
vfoptions0=vfoptionsA;
vfoptions0.e_grid=e_grid; vfoptions0.pi_e=pi_e;
[~,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFnA,Params,DiscountFactorParamNames,[],vfoptions0);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy0,n_d,n_a,n_z,N_j,pi_z,Params,simoptions0);
clear Policy0

%% (X1), (X2), (X3)
for ff=0:1 % fastOLG
    transpathoptions=transpathoptionsbaseline;
    transpathoptions.fastOLG=ff;
    simoptionsA.fastOLG=ff;
    simoptionsB.fastOLG=ff;
    [VPathA,PolicyPathA]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_gridT, pi_zT, DiscountFactorParamNames, ReturnFnA, transpathoptions, vfoptionsA);
    [VPathB,PolicyPathB]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPathB, T, V_final, Policy_final, ParamsB, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_zT, DiscountFactorParamNames, ReturnFnB, transpathoptions, vfoptionsB);
    fprintf('Shocks over the path (experienceassetz), fastOLG=%i: (X1) drifting grids vs the drift as parameters, VPath, this should be zero: %.3e \n',ff,gather(max(abs(VPathA(:)-VPathB(:)))))
    fprintf('Shocks over the path (experienceassetz), fastOLG=%i: (X1) drifting grids vs the drift as parameters, PolicyPath, this should be zero: %.3e \n',ff,gather(max(abs(PolicyPathA(:)-PolicyPathB(:)))))
    AgentDistPathA=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPath, PolicyPathA, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_zT, T,Params, transpathoptions, simoptionsA);
    AgentDistPathB=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPathB, PolicyPathB, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_zT, T,ParamsB, transpathoptions, simoptionsB);
    fprintf('Shocks over the path (experienceassetz), fastOLG=%i: (X2) drifting grids vs the drift as parameters, AgentDistPath, this should be zero: %.3e \n',ff,gather(max(abs(AgentDistPathA(:)-AgentDistPathB(:)))))
    if ff==0
        AggVarsPathA=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluateA, AgentDistPathA, PolicyPathA, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_gridT, transpathoptions, simoptionsA);
        AggVarsPathB=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluateB, AgentDistPathB, PolicyPathB, PricePath, ParamPathB, ParamsB, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptions, simoptionsB);
        fprintf('Shocks over the path (experienceassetz): (X3) drifting grids vs the drift as parameters, AggVars Mean of assets and earnings, this should be zero: %.3e \n',max([max(abs(gather(AggVarsPathA.assets.Mean-AggVarsPathB.assets.Mean))),max(abs(gather(AggVarsPathA.earnings.Mean-AggVarsPathB.earnings.Mean)))]))
        % For scale: the drift must actually move things (the plain model is (B) without the parameters on the path)
        VPathP=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, ParamsB, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_zT, DiscountFactorParamNames, ReturnFnB, transpathoptions, vfoptionsB);
        fprintf('Shocks over the path (experienceassetz): for scale, the drift moves VPath by up to %.3e (must not be zero, or (X1)-(X3) are empty) \n',gather(max(abs(VPathA(:)-VPathP(:)))))
        clear VPathP AggVarsPathA AggVarsPathB
    end
    clear VPathA VPathB PolicyPathA PolicyPathB AgentDistPathA AgentDistPathB
end

%% (X4) one iteration of the GE path
GeneralEqmEqns=struct();
GeneralEqmEqns.earningsgap=@(earnings) earnings-0.5;
transpathoptionsG=transpathoptionsbaseline;
transpathoptionsG.maxiter=1;
transpathoptionsG.GEnewprice=2; % Anderson acceleration
transpathoptionsG.GEnewprice2.howtoupdate={'earningsgap','r',0,0.1};
for ff=0:1
    transpathoptionsG.fastOLG=ff;
    % The agent dist in the GE path is always fastOLG=1 here (its default): the experience asset's a2prime
    % probabilities are a grid interpolation layer, which the GE path does not allow with simoptions.fastOLG=0
    simoptionsA.fastOLG=1;
    simoptionsB.fastOLG=1;
    PricePathGA=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_gridT, pi_zT, ReturnFnA, FnsToEvaluateA, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsG, simoptionsA, vfoptionsA);
    PricePathGB=TransitionPath_Case1_FHorz(PricePath, ParamPathB, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_zT, ReturnFnB, FnsToEvaluateB, GeneralEqmEqns, ParamsB, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsG, simoptionsB, vfoptionsB);
    fprintf('Shocks over the path (experienceassetz), fastOLG=%i: (X4) one iteration of the GE path, drifting grids vs the drift as parameters, this should be zero: %.3e \n',ff,max(abs(PricePathGA.r-PricePathGB.r)))
    fprintf('Shocks over the path (experienceassetz), fastOLG=%i: for scale, that one iteration moves r by up to %.3e (must not be zero, or (X4) is empty) \n',ff,max(abs(PricePathGA.r-PricePath.r)))
end

output=struct();

end
