function output=CoreFHorzTPath_fastOLG2(T,n_d,n_a,n_z,n_e,N_j,d_grid,a_grid,z_grid,pi_z,e_grid,pi_e,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline)
% transpathoptions.fastOLG=2: fastOLG for the local search sweeps, slowOLG for the full sweeps (the first of
% a solve, and every seed or verification sweep). No figures. The model is the one in
% CoreFHorzTPath_localsearch_shooting (doPart(10)): a Cobb-Douglas firm supplies both r and w, and the path
% is a null reform solved from a bumped starting guess, so the solver has to travel.
%
% One solve per row of the configuration table below, done twice with local search, by the SAME algorithm
% with the SAME options: once with fastOLG=1, once with fastOLG=2.
%   (1) the two give the same price path and the same general eqm conditions, exactly: the full sweeps
%       differ only in layout, and slowOLG and fastOLG solve the same problem to the same answer (doParts
%       1-5 check that, sweep by sweep)
%   (2) the DISCRIMINATOR for (1): the fastOLG=2 solve did solve its full sweeps age by age. Run with
%       transpathoptions.verbose=1, each such sweep prints one line; the count must be at least one
%       (iteration 1 is always a full sweep), or (1) compared fastOLG=1 with itself
% The table covers shooting and Anderson, the four shock cases (z, no z, z and e, e without z), with and
% without d, GI and divide-and-conquer (honoured on the full sweeps).
% Then, once:
%   (3) ValueFnOnTransPath_Case1_FHorz: its single pass is a full sweep, so fastOLG=2 must give exactly
%       fastOLG=0; with local search forced (vfoptions.localsearch=2) exactly fastOLG=1
%   (4) the refusals, audited AS A SET: fastOLG=2 without local search; fastOLG=2 with permanent types

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

% columns: fastOLG (not used here), shocks (1=z,2=noz,3=e,4=noz_e), hasd, GI, GEnewprice, divideandconquer
config=[1,1,1,0,3,0;   % shooting, z
        1,1,0,1,3,0;
        1,2,1,1,3,0;   % no z
        1,3,1,0,3,1;   % z and e, divide-and-conquer
        1,4,0,0,3,0;   % e without z
        1,1,1,0,2,0;   % Anderson
        1,1,1,1,2,1];
shockname={'z','noz','e','noz_e'};
algoname={'','Anderson','shooting'};

for cc=1:size(config,1)
    ss=config(cc,2); hasd=config(cc,3); useGI=config(cc,4); GEnewprice=config(cc,5); useDC=config(cc,6);
    hasz=(ss==1 || ss==3);
    hase=(ss==3 || ss==4);

    vfoptions=struct();
    simoptions=struct();
    if useGI
        vfoptions.gridinterplayer=1;
        vfoptions.ngridinterp=5;
        simoptions.gridinterplayer=1;
        simoptions.ngridinterp=5;
    end
    if useDC
        vfoptions.divideandconquer=1;
    end
    if hase
        vfoptions.n_e=n_e; vfoptions.e_grid=e_grid; vfoptions.pi_e=pi_e;
        simoptions.n_e=n_e; simoptions.e_grid=e_grid; simoptions.pi_e=pi_e;
    end
    if hasz
        n_z_cc=n_z; z_grid_cc=z_grid; pi_z_cc=pi_z;
    else
        n_z_cc=0; z_grid_cc=[]; pi_z_cc=[];
    end
    if hasd
        n_d_cc=n_d; d_grid_cc=d_grid;
    else
        n_d_cc=0; d_grid_cc=[];
    end
    % zero assets, mid points for the shocks
    if ss==1
        jequaloneDist=zeros([n_a,n_z],'gpuArray'); jequaloneDist(1,ceil(n_z/2))=1;
    elseif ss==2
        jequaloneDist=zeros([n_a,1],'gpuArray'); jequaloneDist(1)=1;
    elseif ss==3
        jequaloneDist=zeros([n_a,n_z,n_e],'gpuArray'); jequaloneDist(1,ceil(n_z/2),ceil(n_e/2))=1;
    else
        jequaloneDist=zeros([n_a,n_e],'gpuArray'); jequaloneDist(1,ceil(n_e/2))=1;
    end
    % Return fn and aggregates for this model. N is efficiency labour: kappa_j times the shocks (times d)
    FnsToEvaluateGE=struct();
    if hasd
        if ss==1
            ReturnFn=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
            FnsToEvaluateGE.K=@(d,aprime,a,z) a;
            FnsToEvaluateGE.N=@(d,aprime,a,z,kappa_j) kappa_j*z*d;
        elseif ss==2
            ReturnFn=@(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_noe_nosemiz(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
            FnsToEvaluateGE.K=@(d,aprime,a) a;
            FnsToEvaluateGE.N=@(d,aprime,a,kappa_j) kappa_j*d;
        elseif ss==3
            ReturnFn=@(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
            FnsToEvaluateGE.K=@(d,aprime,a,z,e) a;
            FnsToEvaluateGE.N=@(d,aprime,a,z,e,kappa_j) kappa_j*z*e*d;
        else
            ReturnFn=@(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
            FnsToEvaluateGE.K=@(d,aprime,a,e) a;
            FnsToEvaluateGE.N=@(d,aprime,a,e,kappa_j) kappa_j*e*d;
        end
    else
        if ss==1
            ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
            FnsToEvaluateGE.K=@(aprime,a,z) a;
            FnsToEvaluateGE.N=@(aprime,a,z,kappa_j) kappa_j*z;
        elseif ss==2
            ReturnFn=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);
            FnsToEvaluateGE.K=@(aprime,a) a;
            FnsToEvaluateGE.N=@(aprime,a,kappa_j) kappa_j;
        elseif ss==3
            ReturnFn=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
            FnsToEvaluateGE.K=@(aprime,a,z,e) a;
            FnsToEvaluateGE.N=@(aprime,a,z,e,kappa_j) kappa_j*z*e;
        else
            ReturnFn=@(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension);
            FnsToEvaluateGE.K=@(aprime,a,e) a;
            FnsToEvaluateGE.N=@(aprime,a,e,kappa_j) kappa_j*e;
        end
    end
    GeneralEqmEqnsGE=struct();
    GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
    GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

    %% Stationary general eqm, then a null-reform path from a bumped starting guess
    % The stationary solve and V_final use the same vfoptions WITHOUT local search (they are not path solves)
    heteroagentoptionsGE=struct();
    [p_eqm,~]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d_cc,n_a,n_z_cc,N_j,0,pi_z_cc,d_grid_cc,a_grid,z_grid_cc,ReturnFn,FnsToEvaluateGE,GeneralEqmEqnsGE,Params,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptionsGE,simoptions,vfoptions);
    ParamsGE=Params;
    ParamsGE.r=p_eqm.r;
    ParamsGE.w=p_eqm.w;
    [V_finalGE,Policy_finalGE]=ValueFnIter_Case1_FHorz(n_d_cc,n_a,n_z_cc,N_j,d_grid_cc,a_grid,z_grid_cc,pi_z_cc,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initialGE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_finalGE,n_d_cc,n_a,n_z_cc,N_j,pi_z_cc,ParamsGE,simoptions);

    PricePathGE.r=p_eqm.r*ones(1,T);
    PricePathGE.w=p_eqm.w*ones(1,T);
    ParamPathGE.sigma=Params.sigma*ones(1,T);
    Tbump=floor(2*T/3);
    PricePathBumped=PricePathGE;
    PricePathBumped.r(1:Tbump)=PricePathGE.r(1:Tbump).*(1+0.01+0.04*linspace(0,1,Tbump));
    PricePathBumped.w(1:Tbump)=PricePathGE.w(1:Tbump).*(1+0.05-0.04*linspace(0,1,Tbump));

    lab=sprintf('fastOLG=2 %s (%s, d%i, GI%i, DC%i)',algoname{GEnewprice},shockname{ss},hasd,useGI,useDC);

    transpathoptionsS=transpathoptionsbaseline;
    transpathoptionsS.GEnewprice=GEnewprice;
    if GEnewprice==3
        transpathoptionsS.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
    else
        transpathoptionsS.GEnewprice2.howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
    end
    transpathoptionsS.maxiter=250;
    transpathoptionsS.verbose=0;
    transpathoptionsS.toleranceGEcondns=1e-4;
    transpathoptionsS.multiGEcriterion=1;
    simoptions.fastOLG=1;
    vfoptionsL=vfoptions;
    vfoptionsL.localsearch=1;
    vfoptionsL.nlocalsearch=1;

    %% (1) fastOLG=1 against fastOLG=2
    transpathoptions1=transpathoptionsS;
    transpathoptions1.fastOLG=1;
    [PricePath1,GEcondnPath1]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d_cc, n_a, n_z_cc, N_j, d_grid_cc,a_grid,z_grid_cc, pi_z_cc, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions1, simoptions, vfoptionsL);
    transpathoptions2=transpathoptionsS;
    transpathoptions2.fastOLG=2;
    transpathoptions2.verbose=1; % so each full sweep solved age by age says so, see (2)
    out2=evalc('[PricePath2,GEcondnPath2]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d_cc, n_a, n_z_cc, N_j, d_grid_cc,a_grid,z_grid_cc, pi_z_cc, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions2, simoptions, vfoptionsL);');
    fprintf('%s: (1) price path against fastOLG=1, r, this should be zero: %.3e \n',lab,max(abs(PricePath2.r-PricePath1.r)))
    fprintf('%s: (1) price path against fastOLG=1, w, this should be zero: %.3e \n',lab,max(abs(PricePath2.w-PricePath1.w)))
    fprintf('%s: (1) general eqm conditions against fastOLG=1, this should be zero: %.3e \n',lab,max(abs(GEcondnPath2.CapitalMarket-GEcondnPath1.CapitalMarket))+max(abs(GEcondnPath2.LabourMarket-GEcondnPath1.LabourMarket)))
    % anti-vacuity: the solve travelled from the bumped guess
    fprintf('%s: (1) the solve moved the price path off its starting guess, should NOT be zero: %.3e \n',lab,max(abs(PricePath1.r-PricePathBumped.r)))
    %% (2) the discriminator
    nfull=numel(strfind(out2,'a full sweep, solved age by age'));
    fprintf('%s: (2) full sweeps solved age by age (one line each, under verbose), should NOT be zero: %i \n',lab,nfull)
end

%% (3) ValueFnOnTransPath_Case1_FHorz with fastOLG=2
% The with-d, with-z model of the setup, without GI, on a moving path: the value fn alone, no general eqm
ReturnFn=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
vfoptions=struct();
[V_final,Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
PricePathV.r=[linspace(0.03,0.05,floor(T/2)),0.05*ones(1,T-floor(T/2))];
ParamPathV.w=[linspace(1.1,1,floor(T/2)),1*ones(1,T-floor(T/2))];
VPathL=cell(1,3); PolicyPathL=cell(1,3);
for fo=0:2
    transpathoptionsV=transpathoptionsbaseline; transpathoptionsV.fastOLG=fo;
    [VPathL{fo+1},PolicyPathL{fo+1}]=ValueFnOnTransPath_Case1_FHorz(PricePathV, ParamPathV, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsV, vfoptions);
end
fprintf('fastOLG=2 ValueFnOnTransPath: (3) the single pass is a full sweep, against fastOLG=0, this should be zero, V: %.3e \n',max(abs(VPathL{3}(:)-VPathL{1}(:))))
fprintf('fastOLG=2 ValueFnOnTransPath: (3) the single pass is a full sweep, against fastOLG=0, this should be zero, Policy: %.3e \n',max(abs(PolicyPathL{3}(:)-PolicyPathL{1}(:))))
vfoptionsF=vfoptions; vfoptionsF.localsearch=2; vfoptionsF.nlocalsearch=3;
for fo=1:2
    transpathoptionsV=transpathoptionsbaseline; transpathoptionsV.fastOLG=fo;
    [VPathL{fo+1},PolicyPathL{fo+1}]=ValueFnOnTransPath_Case1_FHorz(PricePathV, ParamPathV, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsV, vfoptionsF);
end
fprintf('fastOLG=2 ValueFnOnTransPath: (3) with local search forced, against fastOLG=1, this should be zero, V: %.3e \n',max(abs(VPathL{3}(:)-VPathL{2}(:))))
fprintf('fastOLG=2 ValueFnOnTransPath: (3) with local search forced, against fastOLG=1, this should be zero, Policy: %.3e \n',max(abs(PolicyPathL{3}(:)-PolicyPathL{2}(:))))
clear VPathL PolicyPathL

%% (4) the refusals, as a set
jequaloneDist=zeros([n_a,n_z],'gpuArray'); jequaloneDist(1,ceil(n_z/2))=1;
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_final,n_d,n_a,n_z,N_j,pi_z,Params,struct());
FnsToEvaluate.K=@(d,aprime,a,z) a;
GeneralEqmEqns.CapitalMarket=@(r,K) r-(0.3*((K+0.1)^(-0.7))-0.08);
transpathoptionsR=transpathoptionsbaseline;
transpathoptionsR.fastOLG=2;
transpathoptionsR.GEnewprice=3;
transpathoptionsR.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1};
transpathoptionsR.maxiter=1;
transpathoptionsR.verbose=0;
fired=0; msg='(it ran)';
try
    TransitionPath_Case1_FHorz(PricePathV, ParamPathV, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsR, struct(), vfoptions);
catch ME
    msg=ME.message;
    fired=contains(msg,'fastOLG=2') && contains(msg,'localsearch');
end
fprintf('fastOLG=2: (4) without local search, TransitionPath_Case1_FHorz refuses it by name, this should be one: %i \n',fired)
fprintf('   message: %s \n',msg)
fired=0; msg='(it ran)';
try
    TransitionPath_Case1_FHorz_PType(PricePathV, ParamPathV, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, 2, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, {}, transpathoptionsR, struct(), vfoptions);
catch ME
    msg=ME.message;
    fired=contains(msg,'fastOLG=2') && contains(msg,'permanent types');
end
fprintf('fastOLG=2: (4) with permanent types, TransitionPath_Case1_FHorz_PType refuses it by name, this should be one: %i \n',fired)
fprintf('   message: %s \n',msg)

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
