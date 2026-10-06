function output=CoreFHorzTPathPTypeAlgo_GEptype_islands(T,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,figure_c)
% The transition-path algorithms with permanent types when the general eqm conditions are conditional
% on ptype (transpathoptions.GEptype): two 'islands', each ptype with its own Cobb-Douglas firm and its
% own r and w (so four prices, and four general eqm conditions). The islands do not interact, so the
% PType solve is two solo solves without permanent types side by side.
%   (K1) shooting at a fixed 10 iterations is separable by island, so PType vs solo must be exact
%   (K2) shooting and Anderson both reach toleranceGEcondns (flags), and agree with the solo solves on the
%        eqm path. Anderson mixes the iterates of both islands in one history, so it does not take the
%        solo iterates, and only lands near them.
%   (K3) Anderson with memory=0 is exactly the shooting algorithm, with the four prices
% The model is that of CoreFHorzTPathAlgo_nod_z_noe_nosemiz; the islands differ in kappa_j.

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];

jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

FnsToEvaluateGE.K=@(aprime,a,z) a;
FnsToEvaluateGE.N=@(aprime,a,z,kappa_j) kappa_j*z;

GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

heteroagentoptionsGE=struct();

Names_i={'islandA','islandB'};
N_i=2;
kappa_by_island=[1,1.3];

%% Each island: stationary GE, V_final and initial agent dist, and a bumped starting guess
Tbump=floor(2*T/3);
bumpr=0.01+0.04*linspace(0,1,Tbump);
bumpw=0.05-0.04*linspace(0,1,Tbump);
Params_solo=cell(N_i,1); PricePathBumped_solo=cell(N_i,1); PricePathGE_solo=cell(N_i,1);
V_finalPT=struct(); AgentDist_initialPT=struct();
ParamsPT=Params;
ParamsPT.kappa_j=struct();
for ii=1:N_i
    nm=Names_i{ii};
    Params_ii=Params;
    Params_ii.kappa_j=kappa_by_island(ii)*Params.kappa_j;
    ParamsPT.kappa_j.(nm)=Params_ii.kappa_j;
    [p_eqm_ii,~]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluateGE,GeneralEqmEqnsGE,Params_ii,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptionsGE,simoptions,vfoptions);
    fprintf('Stationary GE (%s): r=%2.8f, w=%2.8f \n',nm,p_eqm_ii.r,p_eqm_ii.w)
    Params_ii.r=p_eqm_ii.r;
    Params_ii.w=p_eqm_ii.w;
    Params_solo{ii}=Params_ii;
    [V_finalPT.(nm),Policy_ii]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initialPT.(nm)=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_ii,n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions);
    PricePathGE_solo{ii}.r=p_eqm_ii.r*ones(1,T);
    PricePathGE_solo{ii}.w=p_eqm_ii.w*ones(1,T);
    PricePathBumped_solo{ii}=PricePathGE_solo{ii};
    PricePathBumped_solo{ii}.r(1:Tbump)=PricePathGE_solo{ii}.r(1:Tbump).*(1+bumpr);
    PricePathBumped_solo{ii}.w(1:Tbump)=PricePathGE_solo{ii}.w(1:Tbump).*(1+bumpw);
end
AgentDist_initialPT.ptweights=[0.5;0.5];
PTypeDistParamNames={'ptweights_algo'};
ParamsPT.ptweights_algo=[0.5;0.5];
% PricePath with prices by ptype, as a structure
PricePathBumpedPT.r.islandA=PricePathBumped_solo{1}.r; PricePathBumpedPT.r.islandB=PricePathBumped_solo{2}.r;
PricePathBumpedPT.w.islandA=PricePathBumped_solo{1}.w; PricePathBumpedPT.w.islandB=PricePathBumped_solo{2}.w;
ParamPathGE.sigma=Params.sigma*ones(1,T); % constant: nothing actually changes, this is a null reform

transpathoptionsGE=transpathoptionsbaseline; % fastOLG=1
transpathoptionsGE.verbose=0;
transpathoptionsGE.GEnewprice=3;
transpathoptionsGE.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; ...
                                            'LabourMarket','w',0,0.1};
transpathoptionsGEPT=transpathoptionsGE;
transpathoptionsGEPT.GEptype={'CapitalMarket','LabourMarket'};

%% (K1) shooting, fixed 10 iterations: separable by island
transpathoptionsK1=transpathoptionsGEPT;
transpathoptionsK1.maxiter=10;
PricePathK1PT=TransitionPath_Case1_FHorz_PType(PricePathBumpedPT, ParamPathGE, T, V_finalPT, AgentDist_initialPT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsPT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsK1, simoptions, vfoptions);
transpathoptionsK1solo=transpathoptionsGE;
transpathoptionsK1solo.maxiter=10;
for ii=1:N_i
    nm=Names_i{ii};
    PricePathK1solo=TransitionPath_Case1_FHorz(PricePathBumped_solo{ii}, ParamPathGE, T, V_finalPT.(nm), AgentDist_initialPT.(nm), jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params_solo{ii}, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsK1solo, simoptions, vfoptions);
    fprintf('GEptype islands, shooting at maxiter=10, PType vs solo (%s), this should be zero, r: %.3e, w: %.3e \n',nm,max(abs(PricePathK1PT.r.(nm)(:)-PricePathK1solo.r(:))),max(abs(PricePathK1PT.w.(nm)(:)-PricePathK1solo.w(:))))
end

%% (K3) Anderson with memory=0 is the shooting algorithm, with prices by ptype
transpathoptionsK3=rmfield(transpathoptionsK1,'GEnewprice3');
transpathoptionsK3.GEnewprice=2;
transpathoptionsK3.GEnewprice2.howtoupdate=transpathoptionsK1.GEnewprice3.howtoupdate;
transpathoptionsK3.anderson.memory=0;
PricePathK3PT=TransitionPath_Case1_FHorz_PType(PricePathBumpedPT, ParamPathGE, T, V_finalPT, AgentDist_initialPT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsPT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsK3, simoptions, vfoptions);
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('GEptype islands, Anderson with memory=0 vs shooting, both at maxiter=10 (%s), this should be zero, r: %.3e, w: %.3e \n',nm,max(abs(PricePathK3PT.r.(nm)(:)-PricePathK1PT.r.(nm)(:))),max(abs(PricePathK3PT.w.(nm)(:)-PricePathK1PT.w.(nm)(:))))
end

%% (K2) to convergence: shooting and Anderson, PType and solo
transpathoptionsK2S=transpathoptionsGEPT;
transpathoptionsK2S.maxiter=250;
tic;
[PricePathK2S,GEcondnPathK2S]=TransitionPath_Case1_FHorz_PType(PricePathBumpedPT, ParamPathGE, T, V_finalPT, AgentDist_initialPT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsPT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsK2S, simoptions, vfoptions);
time_K2S=toc;
transpathoptionsK2A=rmfield(transpathoptionsK2S,'GEnewprice3');
transpathoptionsK2A.GEnewprice=2;
transpathoptionsK2A.GEnewprice2.howtoupdate=transpathoptionsK2S.GEnewprice3.howtoupdate;
transpathoptionsK2A.anderson.memory=5;
tic;
[PricePathK2A,GEcondnPathK2A]=TransitionPath_Case1_FHorz_PType(PricePathBumpedPT, ParamPathGE, T, V_finalPT, AgentDist_initialPT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsPT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsK2A, simoptions, vfoptions);
time_K2A=toc;

% GE condn distance across all four conditions, as the solver measures it
distK2S=zeros(1,T-1); distK2A=zeros(1,T-1);
for ii=1:N_i
    nm=Names_i{ii};
    distK2S=distK2S+GEcondnPathK2S.CapitalMarket.(nm).^2+GEcondnPathK2S.LabourMarket.(nm).^2;
    distK2A=distK2A+GEcondnPathK2A.CapitalMarket.(nm).^2+GEcondnPathK2A.LabourMarket.(nm).^2;
end
distK2S=max(sqrt(distK2S)); distK2A=max(sqrt(distK2A));
fprintf('GEptype islands, shooting reached toleranceGEcondns, should be one: %i (GE condn distance %.3e, toleranceGEcondns 1e-4) \n',distK2S<=1e-4,distK2S)
fprintf('GEptype islands, Anderson (memory=5) reached toleranceGEcondns, should be one: %i (GE condn distance %.3e, toleranceGEcondns 1e-4) \n',distK2A<=1e-4,distK2A)

transpathoptionsK2solo=transpathoptionsGE;
transpathoptionsK2solo.maxiter=250;
for ii=1:N_i
    nm=Names_i{ii};
    PricePathK2solo=TransitionPath_Case1_FHorz(PricePathBumped_solo{ii}, ParamPathGE, T, V_finalPT.(nm), AgentDist_initialPT.(nm), jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params_solo{ii}, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsK2solo, simoptions, vfoptions);
    fprintf('GEptype islands, shooting PType vs solo, same eqm path (%s), should be small, r: %2.10f, w: %2.10f \n',nm,max(abs(PricePathK2S.r.(nm)(:)-PricePathK2solo.r(:))),max(abs(PricePathK2S.w.(nm)(:)-PricePathK2solo.w(:))))
    fprintf('GEptype islands, Anderson PType vs shooting solo, same eqm path (%s), should be small, r: %2.10f, w: %2.10f \n',nm,max(abs(PricePathK2A.r.(nm)(:)-PricePathK2solo.r(:))),max(abs(PricePathK2A.w.(nm)(:)-PricePathK2solo.w(:))))
end
fprintf('Runtime, GEptype islands: shooting %2.2f seconds, Anderson (memory=5) %2.2f seconds \n',time_K2S,time_K2A)

%% A graph of the solutions
figure(figure_c);
subplot(2,1,1); plot(1:T,PricePathBumpedPT.r.islandA,1:T,PricePathK2S.r.islandA,1:T,PricePathK2A.r.islandA,1:T,PricePathBumpedPT.r.islandB,1:T,PricePathK2S.r.islandB,1:T,PricePathK2A.r.islandB)
title('r by island: bumped start, shooting, Anderson')
subplot(2,1,2); plot(1:T,PricePathBumpedPT.w.islandA,1:T,PricePathK2S.w.islandA,1:T,PricePathK2A.w.islandA,1:T,PricePathBumpedPT.w.islandB,1:T,PricePathK2S.w.islandB,1:T,PricePathK2A.w.islandB)
title('w by island: bumped start, shooting, Anderson')
legend('A start','A shooting','A Anderson','B start','B shooting','B Anderson')
% saved by the driver, as in every Core bank

output=1;

end
