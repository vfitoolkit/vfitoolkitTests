function output=CoreInfHorzTPathPTypeAlgo_TwoTypes(T,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,usedecision,figure_c)
% Compares the InfHorz transition-path algorithms with permanent types on one model: GEnewprice=2
% (Anderson acceleration) against GEnewprice=3 (shooting). Both solve the same null-reform general eqm
% path from the same bumped starting guess, and both are timed. usedecision=0 is the model of
% CoreInfHorzTPathAlgo_nod_z_noe_nosemiz, usedecision=1 that of CoreInfHorzTPathAlgo_d_z_noe_nosemiz
% (labour d), each with two permanent types that differ in sigma (given as a structure) and share the
% one Cobb-Douglas firm.
%   (A) each algorithm reaches toleranceGEcondns, printed as a flag, and they agree on the path
%   (B) Anderson with memory=0 is exactly the shooting algorithm
%   (C) N_i=2 identical types against the solver without permanent types, for shooting and for
%       Anderson, each at a fixed number of iterations

vfoptions=struct();
simoptions=struct();

if usedecision==0
    modelname='nod_z_noe';
    n_d=0; d_grid=[];
    ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
    FnsToEvaluateGE.K=@(aprime,a,z) a;
    FnsToEvaluateGE.N=@(aprime,a,z) z;
else
    modelname='d_z_noe';
    ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
    FnsToEvaluateGE.K=@(d,aprime,a,z) a;
    FnsToEvaluateGE.N=@(d,aprime,a,z) z*d;
end

%% General equilibrium: solve the stationary GE first, then run a null-reform transition path
% A Cobb-Douglas firm supplies BOTH prices from its first order conditions:
%     r = firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha)) - firmdelta   (MPK less depreciation)
%     w = (1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha))                (MPL)
% K is aggregate assets and N aggregate efficiency units of labour, both pooled across the two types.
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

heteroagentoptionsGE=struct(); % default fminalgo, and no constraints on r or w

% Two types, sigma differs
Names_i={'typeA','typeB'};
PTypeDistParamNames={'ptweights_algo'};
ParamsPT=Params;
ParamsPT.ptweights_algo=[0.5;0.5];
ParamsPT.sigma=struct();
ParamsPT.sigma.typeA=Params.sigma;
ParamsPT.sigma.typeB=1.5*Params.sigma;

[p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsPT, DiscountFactorParamNames, PTypeDistParamNames, {'r','w'}, heteroagentoptionsGE, simoptions, vfoptions);
fprintf('Stationary GE (%s, two types): r=%2.8f, w=%2.8f \n',modelname,p_eqm.r,p_eqm.w)
fprintf('Stationary GE conditions (these should be close to zero): CapitalMarket=%.3e, LabourMarket=%.3e \n',GEcondns.CapitalMarket,GEcondns.LabourMarket)

% The path needs the EQUILIBRIUM V_final and initial dist, else the null reform is not exact
ParamsGE=ParamsPT;
ParamsGE.r=p_eqm.r;
ParamsGE.w=p_eqm.w;
[V_finalGE,Policy_finalGE]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,vfoptions);
AgentDist_initialGE=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy_finalGE,n_d,n_a,n_z,Names_i,pi_z,ParamsGE,simoptions);

PricePathGE.r=p_eqm.r*ones(1,T);
PricePathGE.w=p_eqm.w*ones(1,T);
ParamPathGE.firmA=Params.firmA*ones(1,T); % constant: nothing actually changes, this is a null reform

transpathoptionsGE=struct();
transpathoptionsGE.verbose=0;
% Both residuals are (price - its firm FOC), so a POSITIVE residual means the price is too high and
% must come down: hence add=0 (subtract) for both.
transpathoptionsGE.GEnewprice=3;
transpathoptionsGE.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; ...
                                            'LabourMarket','w',0,0.1};

%% A bumped starting guess, so the algorithms have to actually travel
Tbump=floor(2*T/3);
bumpr=0.01+0.04*linspace(0,1,Tbump);
bumpw=0.05-0.04*linspace(0,1,Tbump);
PricePathBumped=PricePathGE;
PricePathBumped.r(1:Tbump)=PricePathGE.r(1:Tbump).*(1+bumpr);
PricePathBumped.w(1:Tbump)=PricePathGE.w(1:Tbump).*(1+bumpw);
fprintf('Bumped starting guess, max initial deviation from p_eqm: r %2.8f, w %2.8f \n',max(abs(PricePathBumped.r-p_eqm.r)),max(abs(PricePathBumped.w-p_eqm.w)))

%% (B) GEnewprice=3: the shooting algorithm
transpathoptionsB=transpathoptionsGE;
transpathoptionsB.maxiter=250;
tic;
[PricePathB,GEcondnPathB]=TransitionPath_InfHorz_PType(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsB, simoptions, vfoptions);
time_B=toc;

%% (F) GEnewprice=2: Anderson acceleration of the same shooting map as (B)
transpathoptionsF=rmfield(transpathoptionsB,'GEnewprice3');
transpathoptionsF.GEnewprice=2;
transpathoptionsF.GEnewprice2.howtoupdate=transpathoptionsB.GEnewprice3.howtoupdate;

% (F1) memory=0: every step is a plain shooting step, so Anderson is then exactly the shooting algorithm.
% Both run a fixed 10 iterations, short of convergence.
transpathoptionsF1=transpathoptionsF;
transpathoptionsF1.maxiter=10;
transpathoptionsF1.anderson.memory=0;
PricePathF0=TransitionPath_InfHorz_PType(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsF1, simoptions, vfoptions);
transpathoptionsB1=transpathoptionsB;
transpathoptionsB1.maxiter=10;
PricePathB1=TransitionPath_InfHorz_PType(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsB1, simoptions, vfoptions);

% (F2) And now with a history to work with: safeguard off and on
safeopts=[0,1];
nanderson=length(safeopts);
time_F=zeros(1,nanderson); distF=zeros(1,nanderson);
devF_r=zeros(1,nanderson); devF_w=zeros(1,nanderson);
andersonname=cell(1,nanderson);
for mm=1:nanderson
    andersonname{mm}=sprintf('memory=5, safeguard=%i',safeopts(mm));
    transpathoptionsF.anderson.memory=5;
    transpathoptionsF.anderson.safeguard=safeopts(mm);
    tic;
    [PricePathF,GEcondnPathF]=TransitionPath_InfHorz_PType(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsF, simoptions, vfoptions);
    time_F(mm)=toc;
    distF(mm)=max(sqrt(gather(GEcondnPathF.CapitalMarket(1:T-1)).^2+gather(GEcondnPathF.LabourMarket(1:T-1)).^2));
    devF_r(mm)=max(abs(gather(PricePathF.r)-gather(PricePathB.r)));
    devF_w(mm)=max(abs(gather(PricePathF.w)-gather(PricePathB.w)));
end

%% (A) Each algorithm must reach its own stated tolerance (a flag, so that a run that stops short fails)
distB=max(sqrt(gather(GEcondnPathB.CapitalMarket(1:T-1)).^2+gather(GEcondnPathB.LabourMarket(1:T-1)).^2));
fprintf('Shooting (%s) reached toleranceGEcondns, should be one: %i (GE condn distance %.3e, toleranceGEcondns 1e-4) \n',modelname,distB<=1e-4,distB)
for mm=1:nanderson
    fprintf('Anderson [%s] (%s) reached toleranceGEcondns, should be one: %i (GE condn distance %.3e, toleranceGEcondns 1e-4) \n',andersonname{mm},modelname,distF(mm)<=1e-4,distF(mm))
end
% And they must agree on where the equilibrium path is. Not exactly: each stops the first iteration its
% own residual drops under the tolerance, so they land at different points inside it.
fprintf('Anderson with memory=0 vs shooting (%s), both at maxiter=10, the same algorithm, should be zero, r: %.3e, w: %.3e \n',modelname,max(abs(gather(PricePathF0.r)-gather(PricePathB1.r))),max(abs(gather(PricePathF0.w)-gather(PricePathB1.w))))
for mm=1:nanderson
    fprintf('Anderson [%s] vs shooting (%s), same eqm path, should be small, r: %2.10f, w: %2.10f \n',andersonname{mm},modelname,devF_r(mm),devF_w(mm))
end
fprintf('Runtime (%s), shooting %2.2f seconds \n',modelname,time_B)
for mm=1:nanderson
    fprintf('Runtime (%s), Anderson [%s]: %2.2f seconds \n',modelname,andersonname{mm},time_F(mm))
end

%% (C) N_i=2 identical types against the solver without permanent types
% Both types get the sigma of typeA, weights 0.6 and 0.4, so the pooled aggregates are those of one type,
% and each algorithm must take the same iterates as without permanent types. Stationary eqm and the path
% objects are those of the single type.
[p_eqm1,~]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptions, vfoptions);
Params1=Params;
Params1.r=p_eqm1.r;
Params1.w=p_eqm1.w;
[V_final1,Policy_final1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial1=StationaryDist_InfHorz(Policy_final1,n_d,n_a,n_z,pi_z,simoptions,Params1,[]);
V_final1PT.ptype001=V_final1; V_final1PT.ptype002=V_final1;
AgentDist_initial1PT.ptype001=AgentDist_initial1; AgentDist_initial1PT.ptype002=AgentDist_initial1;
AgentDist_initial1PT.ptweights=[0.6;0.4];
Params1PT=Params1;
Params1PT.ptweights_algo=[0.6;0.4];
PricePath1.r=p_eqm1.r*ones(1,T);
PricePath1.w=p_eqm1.w*ones(1,T);
PricePath1.r(1:Tbump)=PricePath1.r(1:Tbump).*(1+bumpr);
PricePath1.w(1:Tbump)=PricePath1.w(1:Tbump).*(1+bumpw);

PricePathC1_NoPT=TransitionPath_InfHorz(PricePath1, ParamPathGE, T, V_final1, AgentDist_initial1, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params1, DiscountFactorParamNames, transpathoptionsB1, simoptions, vfoptions, []);
PricePathC1_PT=TransitionPath_InfHorz_PType(PricePath1, ParamPathGE, T, V_final1PT, AgentDist_initial1PT, n_d,n_a,n_z, 2, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params1PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsB1, simoptions, vfoptions);
fprintf('Shooting (%s) at maxiter=10, N_i=2 identical types vs no-PType, this should be zero, r: %.3e, w: %.3e \n',modelname,max(abs(gather(PricePathC1_PT.r)-gather(PricePathC1_NoPT.r))),max(abs(gather(PricePathC1_PT.w)-gather(PricePathC1_NoPT.w))))
transpathoptionsC2=transpathoptionsF; % Anderson
transpathoptionsC2.maxiter=10;
transpathoptionsC2.anderson.memory=5;
transpathoptionsC2.anderson.safeguard=0;
PricePathC2_NoPT=TransitionPath_InfHorz(PricePath1, ParamPathGE, T, V_final1, AgentDist_initial1, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params1, DiscountFactorParamNames, transpathoptionsC2, simoptions, vfoptions, []);
PricePathC2_PT=TransitionPath_InfHorz_PType(PricePath1, ParamPathGE, T, V_final1PT, AgentDist_initial1PT, n_d,n_a,n_z, 2, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params1PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsC2, simoptions, vfoptions);
fprintf('Anderson (%s, memory=5) at maxiter=10, N_i=2 identical types vs no-PType, this should be zero, r: %.3e, w: %.3e \n',modelname,max(abs(gather(PricePathC2_PT.r)-gather(PricePathC2_NoPT.r))),max(abs(gather(PricePathC2_PT.w)-gather(PricePathC2_NoPT.w))))
fprintf('Anderson vs shooting (%s) at maxiter=10, no-PType, for scale, this should NOT be zero, r: %.3e \n',modelname,max(abs(gather(PricePathC2_NoPT.r)-gather(PricePathC1_NoPT.r))))

%% Figure: the price paths
figure(figure_c);
clf;
subplot(2,1,1); plot(1:T,gather(PricePathBumped.r),':',1:T,gather(PricePathB.r),'-',1:T,gather(PricePathF.r),'--')
title(['r (',modelname,', two types)']); legend('start','shooting','Anderson')
subplot(2,1,2); plot(1:T,gather(PricePathBumped.w),':',1:T,gather(PricePathB.w),'-',1:T,gather(PricePathF.w),'--')
title(['w (',modelname,', two types)'])

output=struct();

end
