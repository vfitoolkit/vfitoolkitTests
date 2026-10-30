function output=CoreInfHorzTPathPTypeAlgo_GEptype_islands(T,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,figure_c)
% Two 'islands': two permanent types (sigma differs) that each have their own Cobb-Douglas firm, their
% own r and w, and general eqm conditions by ptype (transpathoptions.GEptype). The islands do not
% interact, so the PType solve is two separate problems solved at once.
%   (1) shooting is separable by island, so at a fixed number of iterations it must equal two solo
%       solves without permanent types exactly
%   (2) Anderson with memory=0 is the shooting algorithm
%   (3) shooting and Anderson both converge (a flag each), and to the solo eqm paths
% Each island starts at its own stationary eqm, solved without permanent types, and the path is a null
% reform from a bumped starting guess.

n_d=0;
d_grid=[];
vfoptions=struct();
simoptions=struct();

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
FnsToEvaluateGE.K=@(aprime,a,z) a;
FnsToEvaluateGE.N=@(aprime,a,z) z;

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;
GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));
heteroagentoptionsGE=struct();

Names_i={'islandA','islandB'};
N_i=2;
sigma_i=[Params.sigma,1.5*Params.sigma];
PTypeDistParamNames={'ptweights_algo'};
ParamsPT=Params;
ParamsPT.ptweights_algo=[0.5;0.5];
ParamsPT.sigma=struct();

Tbump=floor(2*T/3);
bumpr=0.01+0.04*linspace(0,1,Tbump);
bumpw=0.05-0.04*linspace(0,1,Tbump);

%% Each island: its stationary eqm (solo), V_final and AgentDist_initial, and its bumped starting guess
V_final_PT=struct(); AgentDist_initial_PT=struct();
PricePath_PT=struct();
PricePath_PT.r=struct(); PricePath_PT.w=struct();
Params_solo=cell(N_i,1); PricePath_solo=cell(N_i,1); p_eqm=cell(N_i,1);
for ii=1:N_i
    nm=Names_i{ii};
    ParamsPT.sigma.(nm)=sigma_i(ii);
    Params_solo{ii}=Params;
    Params_solo{ii}.sigma=sigma_i(ii);
    p_eqm{ii}=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params_solo{ii}, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptions, vfoptions);
    fprintf('Island %s, stationary GE: r=%2.8f, w=%2.8f \n',nm,p_eqm{ii}.r,p_eqm{ii}.w)
    Params_solo{ii}.r=p_eqm{ii}.r;
    Params_solo{ii}.w=p_eqm{ii}.w;
    [V_final_PT.(nm),Policy_ii]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_solo{ii},DiscountFactorParamNames,[],vfoptions);
    AgentDist_initial_PT.(nm)=StationaryDist_InfHorz(Policy_ii,n_d,n_a,n_z,pi_z,simoptions,Params_solo{ii},[]);
    PricePath_solo{ii}.r=p_eqm{ii}.r*ones(1,T);
    PricePath_solo{ii}.w=p_eqm{ii}.w*ones(1,T);
    PricePath_solo{ii}.r(1:Tbump)=PricePath_solo{ii}.r(1:Tbump).*(1+bumpr);
    PricePath_solo{ii}.w(1:Tbump)=PricePath_solo{ii}.w(1:Tbump).*(1+bumpw);
    PricePath_PT.r.(nm)=PricePath_solo{ii}.r;
    PricePath_PT.w.(nm)=PricePath_solo{ii}.w;
end
AgentDist_initial_PT.ptweights=[0.5;0.5];
ParamPathGE.firmA=Params.firmA*ones(1,T); % constant: nothing actually changes, this is a null reform

transpathoptions=struct();
transpathoptions.verbose=0;
transpathoptions.GEnewprice=3;
transpathoptions.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
transpathoptions.GEptype={'CapitalMarket','LabourMarket'};
transpathoptions_solo=rmfield(transpathoptions,'GEptype');

%% (1) Shooting at a fixed number of iterations: PType vs two solo solves
transpathoptions10=transpathoptions;
transpathoptions10.maxiter=10;
transpathoptions_solo10=transpathoptions_solo;
transpathoptions_solo10.maxiter=10;
PricePath_PT10=TransitionPath_InfHorz_PType(PricePath_PT, ParamPathGE, T, V_final_PT, AgentDist_initial_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsPT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions10, simoptions, vfoptions);
for ii=1:N_i
    nm=Names_i{ii};
    PricePath_s10=TransitionPath_InfHorz(PricePath_solo{ii}, ParamPathGE, T, V_final_PT.(nm), AgentDist_initial_PT.(nm), n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params_solo{ii}, DiscountFactorParamNames, transpathoptions_solo10, simoptions, vfoptions, []);
    fprintf('Islands, shooting at maxiter=10, PType vs solo (island %s), this should be zero, r: %.3e, w: %.3e \n',nm,max(abs(gather(PricePath_PT10.r.(nm))-gather(PricePath_s10.r))),max(abs(gather(PricePath_PT10.w.(nm))-gather(PricePath_s10.w))))
end

%% (2) Anderson with memory=0 is the shooting algorithm
transpathoptionsA=rmfield(transpathoptions,'GEnewprice3');
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=transpathoptions.GEnewprice3.howtoupdate;
transpathoptionsA10=transpathoptionsA;
transpathoptionsA10.maxiter=10;
transpathoptionsA10.anderson.memory=0;
PricePath_A10=TransitionPath_InfHorz_PType(PricePath_PT, ParamPathGE, T, V_final_PT, AgentDist_initial_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsPT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsA10, simoptions, vfoptions);
maxdiff=0;
for ii=1:N_i
    nm=Names_i{ii};
    maxdiff=max([maxdiff,max(abs(gather(PricePath_A10.r.(nm))-gather(PricePath_PT10.r.(nm)))),max(abs(gather(PricePath_A10.w.(nm))-gather(PricePath_PT10.w.(nm))))]);
end
fprintf('Islands, Anderson with memory=0 vs shooting, both at maxiter=10, both islands, this should be zero: %.3e \n',maxdiff)

%% (3) Both to convergence, and to the solo eqm paths
transpathoptionsB=transpathoptions;
transpathoptionsB.maxiter=250;
[PricePath_B,GEcondnPath_B]=TransitionPath_InfHorz_PType(PricePath_PT, ParamPathGE, T, V_final_PT, AgentDist_initial_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsPT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsB, simoptions, vfoptions);
transpathoptionsF=transpathoptionsA;
transpathoptionsF.maxiter=250;
transpathoptionsF.anderson.memory=5;
transpathoptionsF.anderson.safeguard=0;
[PricePath_F,GEcondnPath_F]=TransitionPath_InfHorz_PType(PricePath_PT, ParamPathGE, T, V_final_PT, AgentDist_initial_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsPT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsF, simoptions, vfoptions);
transpathoptions_soloB=transpathoptions_solo;
transpathoptions_soloB.maxiter=250;
distB=0; distF=0;
for ii=1:N_i
    nm=Names_i{ii};
    distB=max(distB,max(sqrt(gather(GEcondnPath_B.CapitalMarket.(nm)(1:T-1)).^2+gather(GEcondnPath_B.LabourMarket.(nm)(1:T-1)).^2)));
    distF=max(distF,max(sqrt(gather(GEcondnPath_F.CapitalMarket.(nm)(1:T-1)).^2+gather(GEcondnPath_F.LabourMarket.(nm)(1:T-1)).^2)));
end
fprintf('Islands, shooting reached toleranceGEcondns, should be one: %i (GE condn distance %.3e, toleranceGEcondns 1e-4) \n',distB<=1e-4,distB)
fprintf('Islands, Anderson (memory=5) reached toleranceGEcondns, should be one: %i (GE condn distance %.3e, toleranceGEcondns 1e-4) \n',distF<=1e-4,distF)
for ii=1:N_i
    nm=Names_i{ii};
    PricePath_sB=TransitionPath_InfHorz(PricePath_solo{ii}, ParamPathGE, T, V_final_PT.(nm), AgentDist_initial_PT.(nm), n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params_solo{ii}, DiscountFactorParamNames, transpathoptions_soloB, simoptions, vfoptions, []);
    fprintf('Islands, shooting PType vs solo (island %s), same eqm path, should be small, r: %2.10f, w: %2.10f \n',nm,max(abs(gather(PricePath_B.r.(nm))-gather(PricePath_sB.r))),max(abs(gather(PricePath_B.w.(nm))-gather(PricePath_sB.w))))
    fprintf('Islands, Anderson PType vs solo shooting (island %s), same eqm path, should be small, r: %2.10f, w: %2.10f \n',nm,max(abs(gather(PricePath_F.r.(nm))-gather(PricePath_sB.r))),max(abs(gather(PricePath_F.w.(nm))-gather(PricePath_sB.w))))
end

%% Figure: the interest rate paths of the two islands
figure(figure_c);
clf;
for ii=1:N_i
    nm=Names_i{ii};
    subplot(N_i,1,ii); plot(1:T,gather(PricePath_PT.r.(nm)),':',1:T,gather(PricePath_B.r.(nm)),'-',1:T,gather(PricePath_F.r.(nm)),'--')
    title(['r, ',nm]); legend('start','shooting','Anderson')
end

output=struct();

end
