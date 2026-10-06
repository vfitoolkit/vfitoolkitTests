function output=CoreFHorzTPathPTypeAlgo_nod_z_noe_nosemiz(T,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,figure_c)
% Compares the transition-path algorithms with permanent types on one FHorz model: GEnewprice=2
% (Anderson acceleration) against GEnewprice=3 (shooting), the latter with and without the
% additionalfactor ramp. (Quasi-Newton, GEnewprice=1, is not implemented with permanent types.) All
% of them solve the same null-reform general eqm path from the same bumped starting guess, and all
% are timed. The model is that of CoreFHorzTPathAlgo_nod_z_noe_nosemiz, with two permanent types that
% differ in kappa_j (given as a structure) and share the one Cobb-Douglas firm.
%   (A) each algorithm reaches toleranceGEcondns, printed as a flag, and they agree on the path
%   (B) Anderson with memory=0 is exactly the shooting algorithm
%   (C) N_i=2 identical types against the solver without permanent types, for shooting and for
%       Anderson, each at a fixed number of iterations
%   (D) the order of the howtoupdate rows does not matter
%   (I) fastOLG=0 against fastOLG=1, for shooting and for Anderson

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];

% zero assets, mid point for the shock
jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

%% General equilibrium: solve the stationary GE first, then run a null-reform transition path
% A Cobb-Douglas firm supplies BOTH prices from its first order conditions:
%     r = firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha)) - firmdelta   (MPK less depreciation)
%     w = (1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha))                (MPL)
% K is aggregate assets and N aggregate efficiency units of labour, both pooled across the two types.
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

FnsToEvaluateGE.K=@(aprime,a,z) a;
FnsToEvaluateGE.N=@(aprime,a,z,kappa_j) kappa_j*z;

GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

heteroagentoptionsGE=struct(); % default fminalgo, and no constraints on r or w

% Two types, kappa_j differs
Names_i={'typeA','typeB'};
N_i=2;
PTypeDistParamNames={'ptweights_algo'};
ParamsPT=Params;
ParamsPT.ptweights_algo=[0.5;0.5];
ParamsPT.kappa_j=struct();
ParamsPT.kappa_j.typeA=Params.kappa_j;
ParamsPT.kappa_j.typeB=1.3*Params.kappa_j;

[p_eqm,GEcondns]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, 0, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsPT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, {'r','w'}, heteroagentoptionsGE, simoptions, vfoptions);
fprintf('Stationary GE (two types): r=%2.8f, w=%2.8f \n',p_eqm.r,p_eqm.w)
fprintf('Stationary GE conditions (these should be close to zero): CapitalMarket=%.3e, LabourMarket=%.3e \n',GEcondns.CapitalMarket,GEcondns.LabourMarket)

% The path needs the EQUILIBRIUM V_final and initial dist, else the null reform is not exact
ParamsGE=ParamsPT;
ParamsGE.r=p_eqm.r;
ParamsGE.w=p_eqm.w;
[V_finalGE,Policy_finalGE]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,vfoptions);
AgentDist_initialGE=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_finalGE,n_d,n_a,n_z,N_j,Names_i,pi_z,ParamsGE,simoptions);

PricePathGE.r=p_eqm.r*ones(1,T);
PricePathGE.w=p_eqm.w*ones(1,T);
ParamPathGE.sigma=Params.sigma*ones(1,T); % constant: nothing actually changes, this is a null reform

transpathoptionsGE=transpathoptionsbaseline; % fastOLG=1
transpathoptionsGE.maxiter=25;
transpathoptionsGE.verbose=0;
% Both residuals are (price - its firm FOC), so a POSITIVE residual means the price is too high and
% must come down: hence add=0 (subtract) for both.
transpathoptionsGE.GEnewprice=3;
transpathoptionsGE.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; ...
                                            'LabourMarket','w',0,0.1};

%% A bumped starting guess, so all the algorithms have to actually travel
Tbump=floor(2*T/3);
bumpr=0.01+0.04*linspace(0,1,Tbump);
bumpw=0.05-0.04*linspace(0,1,Tbump);
PricePathBumped=PricePathGE;
PricePathBumped.r(1:Tbump)=PricePathGE.r(1:Tbump).*(1+bumpr);
PricePathBumped.w(1:Tbump)=PricePathGE.w(1:Tbump).*(1+bumpw);
fprintf('Bumped starting guess, max initial deviation from p_eqm: r %2.8f, w %2.8f \n',max(abs(PricePathBumped.r-p_eqm.r)),max(abs(PricePathBumped.w-p_eqm.w)))

%% (B) GEnewprice=3: the shooting algorithm, no ramp
transpathoptionsB=transpathoptionsGE;
transpathoptionsB.maxiter=250;
tic;
[PricePathB,GEcondnPathB]=TransitionPath_Case1_FHorz_PType(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsB, simoptions, vfoptions);
time_B=toc;

%% (C) GEnewprice=3 again, this time with the additionalfactor ramp
transpathoptionsC=transpathoptionsB;
transpathoptionsC.GEnewprice3.additionalfactor=[3,10,30];
tic;
[PricePathC,GEcondnPathC]=TransitionPath_Case1_FHorz_PType(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsC, simoptions, vfoptions);
time_C=toc;

%% (F) GEnewprice=2: Anderson acceleration of the same shooting map as (B)
transpathoptionsF=rmfield(transpathoptionsB,'GEnewprice3');
transpathoptionsF.GEnewprice=2;
transpathoptionsF.GEnewprice2.howtoupdate=transpathoptionsB.GEnewprice3.howtoupdate;

% (F1) memory=0: every step is a plain shooting step, so Anderson is then exactly the shooting algorithm.
% Both run a fixed 10 iterations, short of convergence.
transpathoptionsF1=transpathoptionsF;
transpathoptionsF1.maxiter=10;
transpathoptionsF1.anderson.memory=0;
[PricePathF0,~]=TransitionPath_Case1_FHorz_PType(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsF1, simoptions, vfoptions);
transpathoptionsB1=transpathoptionsB;
transpathoptionsB1.maxiter=10;
[PricePathB1,~]=TransitionPath_Case1_FHorz_PType(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsB1, simoptions, vfoptions);

% (F2) And now with a history to work with. The same five configurations as the bank without permanent types.
memopts=[3,5,10,5,5];
safeopts=[0,0,0,1,0];
typeopts={'II','II','II','II','I'};
nanderson=length(memopts);
time_F=zeros(1,nanderson); distF=zeros(1,nanderson);
devF_r=zeros(1,nanderson); devF_w=zeros(1,nanderson);
andersonname=cell(1,nanderson);
for mm=1:nanderson
    andersonname{mm}=sprintf('memory=%2i, safeguard=%i, type=%s',memopts(mm),safeopts(mm),typeopts{mm});
    transpathoptionsF.anderson.memory=memopts(mm);
    transpathoptionsF.anderson.safeguard=safeopts(mm);
    transpathoptionsF.anderson.type=typeopts{mm};
    fprintf('--- Anderson acceleration, %s --- \n',andersonname{mm})
    tic;
    [PricePathF,GEcondnPathF]=TransitionPath_Case1_FHorz_PType(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsF, simoptions, vfoptions);
    time_F(mm)=toc;
    distF(mm)=max(sqrt(GEcondnPathF.CapitalMarket.^2+GEcondnPathF.LabourMarket.^2));
    devF_r(mm)=max(abs(PricePathF.r-PricePathB.r));
    devF_w(mm)=max(abs(PricePathF.w-PricePathB.w));
end

%% (A) Each algorithm must reach its own stated tolerance (a flag, so that a run that stops short fails)
distB=max(sqrt(GEcondnPathB.CapitalMarket.^2+GEcondnPathB.LabourMarket.^2));
distC=max(sqrt(GEcondnPathC.CapitalMarket.^2+GEcondnPathC.LabourMarket.^2));
fprintf('Shooting reached toleranceGEcondns, should be one: %i (GE condn distance %.3e, toleranceGEcondns 1e-4) \n',distB<=1e-4,distB)
fprintf('Shooting+ramp reached toleranceGEcondns, should be one: %i (GE condn distance %.3e, toleranceGEcondns 1e-4) \n',distC<=1e-4,distC)
for mm=1:nanderson
    fprintf('Anderson [%s] reached toleranceGEcondns, should be one: %i (GE condn distance %.3e, toleranceGEcondns 1e-4) \n',andersonname{mm},distF(mm)<=1e-4,distF(mm))
end

%% And they must agree on where the equilibrium path is. Not exactly: each stops the first iteration its
% own residual drops under the tolerance, so they land at different points inside it.
fprintf('Shooting with/without ramp, same eqm path, should be small, r: %2.10f, w: %2.10f \n',max(abs(PricePathC.r-PricePathB.r)),max(abs(PricePathC.w-PricePathB.w)))
fprintf('Anderson with memory=0 vs shooting, both at maxiter=10, the same algorithm, should be zero, r: %.3e, w: %.3e \n',max(abs(PricePathF0.r-PricePathB1.r)),max(abs(PricePathF0.w-PricePathB1.w)))
for mm=1:nanderson
    fprintf('Anderson [%s] vs shooting, same eqm path, should be small, r: %2.10f, w: %2.10f \n',andersonname{mm},devF_r(mm),devF_w(mm))
end
fprintf('Runtime, shooting %2.2f seconds, shooting+ramp %2.2f seconds \n',time_B,time_C)
for mm=1:nanderson
    fprintf('Runtime, Anderson [%s]: %2.2f seconds \n',andersonname{mm},time_F(mm))
end

%% (C) N_i=2 identical types against the solver without permanent types
% Both types get the kappa_j of typeA, weights 0.6 and 0.4, so the pooled aggregates are those of one
% type, and each algorithm must take the same iterates as without permanent types. Stationary eqm and
% the path objects are those of the single type.
[p_eqm1,~]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluateGE,GeneralEqmEqnsGE,Params,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptionsGE,simoptions,vfoptions);
Params1=Params;
Params1.r=p_eqm1.r;
Params1.w=p_eqm1.w;
[V_final1,Policy_final1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_final1,n_d,n_a,n_z,N_j,pi_z,Params1,simoptions);
V_final1PT.ptype001=V_final1; V_final1PT.ptype002=V_final1;
AgentDist_initial1PT.ptype001=AgentDist_initial1; AgentDist_initial1PT.ptype002=AgentDist_initial1;
AgentDist_initial1PT.ptweights=[0.6;0.4];
Params1PT=Params1;
Params1PT.ptweights_algo=[0.6;0.4];
PricePath1.r=p_eqm1.r*ones(1,T);
PricePath1.w=p_eqm1.w*ones(1,T);
PricePath1.r(1:Tbump)=PricePath1.r(1:Tbump).*(1+bumpr);
PricePath1.w(1:Tbump)=PricePath1.w(1:Tbump).*(1+bumpw);

transpathoptionsC1=transpathoptionsB1; % shooting, 10 iterations
PricePathC1_NoPT=TransitionPath_Case1_FHorz(PricePath1, ParamPathGE, T, V_final1, AgentDist_initial1, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params1, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsC1, simoptions, vfoptions);
PricePathC1_PT=TransitionPath_Case1_FHorz_PType(PricePath1, ParamPathGE, T, V_final1PT, AgentDist_initial1PT, jequaloneDist, n_d, n_a, n_z, N_j, 2, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params1PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsC1, simoptions, vfoptions);
fprintf('Shooting at maxiter=10, N_i=2 identical types vs no-PType, this should be zero, r: %.3e, w: %.3e \n',max(abs(PricePathC1_PT.r-PricePathC1_NoPT.r)),max(abs(PricePathC1_PT.w-PricePathC1_NoPT.w)))
transpathoptionsC2=transpathoptionsF; % Anderson
transpathoptionsC2.maxiter=10;
transpathoptionsC2.anderson.memory=5;
transpathoptionsC2.anderson.safeguard=0;
transpathoptionsC2.anderson.type='II';
PricePathC2_NoPT=TransitionPath_Case1_FHorz(PricePath1, ParamPathGE, T, V_final1, AgentDist_initial1, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params1, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsC2, simoptions, vfoptions);
PricePathC2_PT=TransitionPath_Case1_FHorz_PType(PricePath1, ParamPathGE, T, V_final1PT, AgentDist_initial1PT, jequaloneDist, n_d, n_a, n_z, N_j, 2, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params1PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsC2, simoptions, vfoptions);
fprintf('Anderson (memory=5) at maxiter=10, N_i=2 identical types vs no-PType, this should be zero, r: %.3e, w: %.3e \n',max(abs(PricePathC2_PT.r-PricePathC2_NoPT.r)),max(abs(PricePathC2_PT.w-PricePathC2_NoPT.w)))

%% (D) The order of the howtoupdate rows must not matter
% Two iterations, as shooting returns the path it evaluated: after one it is the starting guess whatever
% the rows say, so the second path is the first update.
transpathoptionsD=transpathoptionsGE;
transpathoptionsD.maxiter=2;
PricePathD1=TransitionPath_Case1_FHorz_PType(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsD, simoptions, vfoptions);
transpathoptionsD.GEnewprice3.howtoupdate={'LabourMarket','w',0,0.1; ...
                                          'CapitalMarket','r',0,0.1}; % the other order
PricePathD2=TransitionPath_Case1_FHorz_PType(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsD, simoptions, vfoptions);
fprintf('howtoupdate rows in price order vs eqn order, two iterations, this should be zero, r: %.3e, w: %.3e \n',max(abs(PricePathD2.r-PricePathD1.r)),max(abs(PricePathD2.w-PricePathD1.w)))
fprintf('howtoupdate rows, for scale the update moved r by %.3e (must not be zero, or the check above is empty) \n',max(abs(PricePathD1.r-PricePathBumped.r)))

%% (I) fastOLG=0 vs fastOLG=1
transpathoptionsI=transpathoptionsB1; % 10 iterations of shooting
transpathoptionsI.fastOLG=0;
PricePathI0=TransitionPath_Case1_FHorz_PType(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsI, simoptions, vfoptions);
fprintf('Shooting at maxiter=10, fastOLG=0 vs fastOLG=1, this should be zero, r: %.3e, w: %.3e \n',max(abs(PricePathI0.r-PricePathB1.r)),max(abs(PricePathI0.w-PricePathB1.w)))
transpathoptionsI=transpathoptionsF; % Anderson, last configuration of (F2) overwritten below
transpathoptionsI.maxiter=10;
transpathoptionsI.anderson.memory=5;
transpathoptionsI.anderson.safeguard=0;
transpathoptionsI.anderson.type='II';
transpathoptionsI.fastOLG=1;
PricePathIA1=TransitionPath_Case1_FHorz_PType(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsI, simoptions, vfoptions);
transpathoptionsI.fastOLG=0;
PricePathIA0=TransitionPath_Case1_FHorz_PType(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsI, simoptions, vfoptions);
fprintf('Anderson (memory=5) at maxiter=10, fastOLG=0 vs fastOLG=1, this should be zero, r: %.3e, w: %.3e \n',max(abs(PricePathIA0.r-PricePathIA1.r)),max(abs(PricePathIA0.w-PricePathIA1.w)))

%% A graph of the solutions, so that a check that fails can be looked at
figure(figure_c);
subplot(2,1,1); plot(1:T,PricePathBumped.r,1:T,PricePathB.r,1:T,PricePathF.r,1:T,PricePathGE.r,'--')
title('r: bumped start, shooting, last Anderson, eqm')
subplot(2,1,2); plot(1:T,PricePathBumped.w,1:T,PricePathB.w,1:T,PricePathF.w,1:T,PricePathGE.w,'--')
title('w: bumped start, shooting, last Anderson, eqm')
legend('start','shooting','Anderson','eqm')
% saved by the driver, as in every Core bank

output=1;

end
