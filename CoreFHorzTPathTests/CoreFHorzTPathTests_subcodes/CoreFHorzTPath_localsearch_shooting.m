function output=CoreFHorzTPath_localsearch_shooting(T,n_d,n_a,n_z,n_e,N_j,d_grid,a_grid,z_grid,pi_z,e_grid,pi_e,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline)
% vfoptions.localsearch=1 inside a real general eqm transition path. No figures.
% The FHorz counterpart of CoreInfHorzTPathTests doPart(6). The model is the one in
% CoreFHorzTPathAlgoTests: a Cobb-Douglas firm supplies both r and w, and the path is a null reform
% solved from a bumped starting guess, so the solver has to travel.
%
% One solve per row of the configuration table below: each is solved twice, once plainly and once with
% local search, by the SAME algorithm on the SAME layout with the SAME divideandconquer setting:
%   (2) the two must converge to the same price path
%   (3) the local search answer is CERTIFIED: its returned general eqm conditions are within tolerance,
%       and an INDEPENDENT plain sweep at its prices gives the same conditions
% The table covers both layouts, shooting and Anderson, divide-and-conquer (honoured on the seed sweeps),
% and the four shock cases (z, no z, z and e, e without z), each with or without d and GI.
% Then, once, on the with-d slowOLG model without GI:
%   (4) ValueFnOnTransPath_Case1_FHorz ignores localsearch=1 with a warning, and localsearch=2 forces it
%   (5) the refusals, audited AS A SET: each must fire, and with its own message

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

% columns: fastOLG, shocks (1=z,2=noz,3=e,4=noz_e), hasd, GI, GEnewprice, divideandconquer
config=[0,1,0,0,3,0;   % the four slowOLG z corners by shooting (the wave 1 runs)
        0,1,0,1,3,0;
        0,1,1,0,3,0;
        0,1,1,1,3,0;
        1,1,0,0,3,0;   % the four fastOLG z corners by shooting
        1,1,0,1,3,0;
        1,1,1,0,3,0;
        1,1,1,1,3,0;
        0,1,0,0,2,0;   % Anderson
        0,1,1,1,2,0;
        1,1,1,0,2,0;
        0,1,1,0,3,1;   % divide-and-conquer, honoured on the seed sweeps
        1,1,0,1,3,1;
        0,3,1,0,3,0;   % z and e
        1,3,0,1,3,0;
        0,2,0,0,3,0;   % no z
        1,2,1,1,3,0;
        0,4,1,1,3,0;   % e without z
        1,4,0,0,3,0];
layoutname={'slowOLG','fastOLG'};
shockname={'z','noz','e','noz_e'};
algoname={'','Anderson','shooting'};

for cc=1:size(config,1)
    fastOLG=config(cc,1); ss=config(cc,2); hasd=config(cc,3); useGI=config(cc,4); GEnewprice=config(cc,5); useDC=config(cc,6);
    hasz=(ss==1 || ss==3);
    hase=(ss==3 || ss==4);
    lab=sprintf('localsearch FHorz %s (%s, %s, d%i, GI%i, DC%i)',algoname{GEnewprice},layoutname{fastOLG+1},shockname{ss},hasd,useGI,useDC);

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

    transpathoptionsS=transpathoptionsbaseline;
    transpathoptionsS.fastOLG=fastOLG;
    transpathoptionsS.GEnewprice=GEnewprice;
    if GEnewprice==3
        transpathoptionsS.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
    else
        transpathoptionsS.GEnewprice2.howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
    end
    transpathoptionsS.maxiter=250;
    transpathoptionsS.verbose=0;
    % Set here, not left to the default, because the checks in (3) read them: the distance below is the
    % one the solver uses under multiGEcriterion=1 (root sum of squares per period, then max over periods)
    transpathoptionsS.toleranceGEcondns=1e-4;
    transpathoptionsS.multiGEcriterion=1;
    if fastOLG==1
        simoptions.fastOLG=1;
    end

    %% (1) the plain solve
    tic;
    [PricePathB,GEcondnPathB]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d_cc, n_a, n_z_cc, N_j, d_grid_cc,a_grid,z_grid_cc, pi_z_cc, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptions);
    timeB=toc;

    %% (2) the same solve with local search
    vfoptionsL=vfoptions;
    vfoptionsL.localsearch=1;
    vfoptionsL.nlocalsearch=1;
    tic;
    [PricePathL,GEcondnPathL]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d_cc, n_a, n_z_cc, N_j, d_grid_cc,a_grid,z_grid_cc, pi_z_cc, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptionsL);
    timeL=toc;
    fprintf('%s: plain took %.2f seconds, with local search %.2f seconds \n',lab,timeB,timeL)
    % Both are solved to toleranceGEcondns from the same start, so they agree to about that tolerance
    fprintf('%s: (2) price path with local search against without, r, should be close to zero: %.3e \n',lab,max(abs(PricePathL.r-PricePathB.r)))
    fprintf('%s: (2) price path with local search against without, w, should be close to zero: %.3e \n',lab,max(abs(PricePathL.w-PricePathB.w)))

    %% (3) CERTIFIED
    distB=max(sqrt(GEcondnPathB.CapitalMarket.^2+GEcondnPathB.LabourMarket.^2));
    distL=max(sqrt(GEcondnPathL.CapitalMarket.^2+GEcondnPathL.LabourMarket.^2));
    fprintf('%s: (3) plain general eqm conditions within tolerance, should be one: %i \n',lab,distB<=transpathoptionsS.toleranceGEcondns)
    fprintf('%s: (3) local search general eqm conditions within tolerance, should be one: %i \n',lab,distL<=transpathoptionsS.toleranceGEcondns)
    % Independently of the solver's bookkeeping: one plain shooting iteration AT the local search answer
    transpathoptionsX=transpathoptionsS;
    transpathoptionsX.maxiter=1;
    transpathoptionsX.GEnewprice=3;
    transpathoptionsX.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
    [~,GEcondnPathX]=TransitionPath_Case1_FHorz(PricePathL, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d_cc, n_a, n_z_cc, N_j, d_grid_cc,a_grid,z_grid_cc, pi_z_cc, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsX, simoptions, vfoptions);
    distX=max(sqrt(GEcondnPathX.CapitalMarket.^2+GEcondnPathX.LabourMarket.^2));
    fprintf('%s: (3) the local search answer re-evaluated by an independent plain sweep, general eqm conditions within tolerance, should be one: %i \n',lab,distX<=transpathoptionsS.toleranceGEcondns)
    if GEnewprice==3
        % Shooting returns the conditions of its final (standard) sweep at the returned prices, so these
        % agree exactly. Anderson returns the certifying sweep's, which is the same evaluation.
        fprintf('%s: (3) local search conditions against the independent plain sweep at the same prices, this should be zero: %.3e \n',lab,max(abs(GEcondnPathL.CapitalMarket-GEcondnPathX.CapitalMarket))+max(abs(GEcondnPathL.LabourMarket-GEcondnPathX.LabourMarket)))
    end
end

%% (4) and (5) once, on the with-d slowOLG z model without GI
ReturnFn=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
FnsToEvaluateGE=struct();
FnsToEvaluateGE.K=@(d,aprime,a,z) a;
FnsToEvaluateGE.N=@(d,aprime,a,z,kappa_j) kappa_j*z*d;
GeneralEqmEqnsGE=struct();
GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));
jequaloneDist=zeros([n_a,n_z],'gpuArray'); jequaloneDist(1,ceil(n_z/2))=1;
vfoptions=struct();
simoptions=struct();
[p_eqm,~]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluateGE,GeneralEqmEqnsGE,Params,DiscountFactorParamNames,[],[],[],{'r','w'},struct(),simoptions,vfoptions);
ParamsGE=Params;
ParamsGE.r=p_eqm.r;
ParamsGE.w=p_eqm.w;
[V_finalGE,Policy_finalGE]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions);
AgentDist_initialGE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_finalGE,n_d,n_a,n_z,N_j,pi_z,ParamsGE,simoptions);
PricePathGE.r=p_eqm.r*ones(1,T);
PricePathGE.w=p_eqm.w*ones(1,T);
ParamPathGE.sigma=Params.sigma*ones(1,T);
Tbump=floor(2*T/3);
PricePathBumped=PricePathGE;
PricePathBumped.r(1:Tbump)=PricePathGE.r(1:Tbump).*(1+0.01+0.04*linspace(0,1,Tbump));
PricePathBumped.w(1:Tbump)=PricePathGE.w(1:Tbump).*(1+0.05-0.04*linspace(0,1,Tbump));
transpathoptionsS=transpathoptionsbaseline;
transpathoptionsS.fastOLG=0;
transpathoptionsS.GEnewprice=3;
transpathoptionsS.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
transpathoptionsS.verbose=0;
transpathoptionsVF=transpathoptionsbaseline;
transpathoptionsVF.fastOLG=0;
vfoptionsL=struct();
vfoptionsL.localsearch=1;
vfoptionsL.nlocalsearch=1;
lab='localsearch FHorz shooting';

%% (4) ValueFnOnTransPath_Case1_FHorz ignores localsearch=1, with a warning, and localsearch=2 forces it
lastwarn('');
[VPath0,PolicyPath0]=ValueFnOnTransPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, Policy_finalGE, ParamsGE, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsVF, vfoptions);
[VPath1,PolicyPath1]=ValueFnOnTransPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, Policy_finalGE, ParamsGE, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsVF, vfoptionsL);
[warnmsg,~]=lastwarn;
fprintf('%s: (4) ValueFnOnTransPath with localsearch=1 warns that it is ignored, should be one: %i \n',lab,contains(warnmsg,'localsearch=1 is IGNORED'))
fprintf('%s: (4) ValueFnOnTransPath with localsearch=1 against without, V, this should be zero: %.3e \n',lab,max(abs(VPath1(:)-VPath0(:))))
fprintf('%s: (4) ValueFnOnTransPath with localsearch=1 against without, Policy, this should be zero: %.3e \n',lab,max(abs(PolicyPath1(:)-PolicyPath0(:))))
% localsearch=2: forced, around the dispatchers' stay-put reference. The window binds (the young save,
% the old dissave), so V moves, and without GI it can only fall. At the cap it is the plain search again.
vfoptionsF=vfoptionsL; vfoptionsF.localsearch=2;
[VPath2,~]=ValueFnOnTransPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, Policy_finalGE, ParamsGE, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsVF, vfoptionsF);
fprintf('%s: (4) ValueFnOnTransPath with localsearch=2, nlocalsearch=1, V against without, should NOT be zero: %.3e \n',lab,max(abs(VPath2(:)-VPath0(:))))
fprintf('%s: (4) ValueFnOnTransPath with localsearch=2, nlocalsearch=1, V is never above without, this should be zero: %.3e \n',lab,max(max(VPath2(:)-VPath0(:)),0))
vfoptionsF.nlocalsearch=(n_a-1)/2;
[VPath2c,PolicyPath2c]=ValueFnOnTransPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, Policy_finalGE, ParamsGE, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsVF, vfoptionsF);
fprintf('%s: (4) ValueFnOnTransPath with localsearch=2 at the cap, V against without, this should be zero: %.3e \n',lab,max(abs(VPath2c(:)-VPath0(:))))
fprintf('%s: (4) ValueFnOnTransPath with localsearch=2 at the cap, Policy entries that differ (ties only): %i of %i \n',lab,sum(PolicyPath2c(:)~=PolicyPath0(:)),numel(PolicyPath0))
clear VPath0 VPath1 VPath2 VPath2c PolicyPath0 PolicyPath1 PolicyPath2c
err4=0;
vfoptionsE=vfoptionsL; vfoptionsE.localsearch=3;
try
    ValueFnOnTransPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, Policy_finalGE, ParamsGE, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsVF, vfoptionsE);
catch ME
    err4=contains(ME.message,'must be 0, 1 (ignored here, with a warning) or 2 (forced)');
end
fprintf('%s: (4) ValueFnOnTransPath with localsearch=3 is refused by name, should be one: %i \n',lab,err4)
err4=0;
vfoptionsE=vfoptionsL; vfoptionsE.localsearch=2; vfoptionsE.experienceasset=1;
try
    ValueFnOnTransPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, Policy_finalGE, ParamsGE, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsVF, vfoptionsE);
catch ME
    err4=contains(ME.message,'vfoptions.localsearch is not implemented with an experience asset');
end
fprintf('%s: (4) ValueFnOnTransPath with localsearch=2 and an experience asset is refused by name, should be one: %i \n',lab,err4)

%% (5) The refusals in TransitionPath_Case1_FHorz, as a set
% Every one of these runs BEFORE any solving. maxiter=1 all the same, so that a refusal that has been
% lifted costs one iteration rather than a full solve before it shows up as a 0.
transpathoptionsS.maxiter=1;
refusalname={}; refusalfired=[];
% localsearch=2 reaching the TPath command
vfoptionsE=vfoptionsL; vfoptionsE.localsearch=2; fired=0;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptionsE);
catch ME
    fired=contains(ME.message,'must be 0 or 1 here');
end
refusalname{end+1}='localsearch=2'; refusalfired(end+1)=fired;
% nlocalsearch not set
vfoptionsE=vfoptionsL; vfoptionsE=rmfield(vfoptionsE,'nlocalsearch'); fired=0;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptionsE);
catch ME
    fired=contains(ME.message,'you must set vfoptions.nlocalsearch');
end
refusalname{end+1}='nlocalsearch not set'; refusalfired(end+1)=fired;
% a ratchet step below one
vfoptionsE=vfoptionsL; vfoptionsE.nlocalsearchup=0; fired=0;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptionsE);
catch ME
    fired=contains(ME.message,'must both be at least 1');
end
refusalname{end+1}='nlocalsearchup=0'; refusalfired(end+1)=fired;
% a window wider than the grid
vfoptionsE=vfoptionsL; vfoptionsE.nlocalsearch=(n_a+1)/2; fired=0;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptionsE);
catch ME
    fired=contains(ME.message,'too large for the grid');
end
refusalname{end+1}='nlocalsearch too large'; refusalfired(end+1)=fired;
% quasi-Newton
transpathoptionsE=transpathoptionsS; transpathoptionsE.GEnewprice=1; fired=0;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptionsL);
catch ME
    fired=contains(ME.message,'requires transpathoptions.GEnewprice=3 (shooting) or 2 (Anderson)');
end
refusalname{end+1}='GEnewprice=1'; refusalfired(end+1)=fired;
% Anderson type I (refused inside AndersonAccelerationTPath)
transpathoptionsE=rmfield(transpathoptionsS,'GEnewprice3'); transpathoptionsE.GEnewprice=2; transpathoptionsE.GEnewprice2.howtoupdate=transpathoptionsS.GEnewprice3.howtoupdate; transpathoptionsE.anderson.type='I'; fired=0;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptionsL);
catch ME
    fired=contains(ME.message,'anderson.type=''I'' is not implemented with vfoptions.localsearch');
end
refusalname{end+1}='anderson.type=I'; refusalfired(end+1)=fired;
% a semi-exogenous state
vfoptionsE=vfoptionsL; vfoptionsE.n_semiz=2; fired=0;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptionsE);
catch ME
    fired=contains(ME.message,'not implemented with a semi-exogenous state');
end
refusalname{end+1}='semiz'; refusalfired(end+1)=fired;
% an experience asset (Step1's ExpAsset branch bypasses the dispatcher)
vfoptionsE=vfoptionsL; vfoptionsE.experienceasset=1; fired=0;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptionsE);
catch ME
    fired=contains(ME.message,'not implemented with an experience asset');
end
refusalname{end+1}='experienceasset'; refusalfired(end+1)=fired;
% quasi-hyperbolic discounting (Step1's QH branch bypasses the dispatcher)
vfoptionsE=vfoptionsL; vfoptionsE.exoticpreferences='QuasiHyperbolic'; vfoptionsE.QHadditionaldiscount='beta0'; ParamsE=ParamsGE; ParamsE.beta0=1; fired=0;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptionsE);
catch ME
    fired=contains(ME.message,'not implemented with exoticpreferences');
end
refusalname{end+1}='QuasiHyperbolic'; refusalfired(end+1)=fired;
for ee=1:length(refusalname)
    fprintf('%s: (5) local search with %s is refused by name, should be one: %i \n',lab,refusalname{ee},refusalfired(ee))
end

output=1;

end
