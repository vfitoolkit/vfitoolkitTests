function output=CoreFHorzTPath_localsearch_shooting(T,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline)
% vfoptions.localsearch=1 inside a real general eqm transition path, by shooting (GEnewprice=3), finite
% horizon, slowOLG (transpathoptions.fastOLG=0). No figures.
% The FHorz counterpart of CoreInfHorzTPathTests doPart(6). The model is the one in
% CoreFHorzTPathAlgoTests: a Cobb-Douglas firm supplies both r and w, and the path is a null reform
% solved from a bumped starting guess, so the solver has to travel.
%
% For each of the four corners (without/with d, without/with the grid interpolation layer):
%   (1) plain shooting, the baseline
%   (2) shooting with local search at nlocalsearch=1, which must converge to the same price path
%   (3) its answer is CERTIFIED: the returned general eqm conditions come from the final standard sweep
% Then, once, on the with-d model without GI:
%   (4) ValueFnOnTransPath_Case1_FHorz ignores localsearch=1 with a warning, and gives the same answer
%   (5) the refusals, audited AS A SET: each must fire, and with its own message. A refusal that has
%       been lifted, or that some other error now pre-empts, shows up as a 0 here rather than vanishing.

transpathoptionsbaseline.fastOLG=0; % local search is implemented for slowOLG only

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

for cc=1:4
    hasd=(cc>=3);
    useGI=(cc==2 || cc==4);
    if hasd
        n_d_cc=n_d;
        d_grid_cc=d_grid;
        ReturnFn=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
        FnsToEvaluateGE=struct();
        FnsToEvaluateGE.K=@(d,aprime,a,z) a;
        FnsToEvaluateGE.N=@(d,aprime,a,z,kappa_j) kappa_j*z*d;
    else
        n_d_cc=0;
        d_grid_cc=[];
        ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
        FnsToEvaluateGE=struct();
        FnsToEvaluateGE.K=@(aprime,a,z) a;
        FnsToEvaluateGE.N=@(aprime,a,z,kappa_j) kappa_j*z;
    end
    GeneralEqmEqnsGE=struct();
    GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
    GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));
    if hasd
        lab=sprintf('localsearch FHorz shooting (d, GI%i)',useGI);
    else
        lab=sprintf('localsearch FHorz shooting (nod, GI%i)',useGI);
    end

    vfoptions=struct();
    simoptions=struct();
    if useGI
        vfoptions.gridinterplayer=1;
        vfoptions.ngridinterp=5;
        simoptions.gridinterplayer=1;
        simoptions.ngridinterp=5;
    end

    jequaloneDist=zeros([n_a,n_z],'gpuArray');
    jequaloneDist(1,ceil(n_z/2))=1;

    %% Stationary general eqm, then a null-reform path from a bumped starting guess
    heteroagentoptionsGE=struct();
    [p_eqm,~]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d_cc,n_a,n_z,N_j,0,pi_z,d_grid_cc,a_grid,z_grid,ReturnFn,FnsToEvaluateGE,GeneralEqmEqnsGE,Params,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptionsGE,simoptions,vfoptions);
    ParamsGE=Params;
    ParamsGE.r=p_eqm.r;
    ParamsGE.w=p_eqm.w;
    [V_finalGE,Policy_finalGE]=ValueFnIter_Case1_FHorz(n_d_cc,n_a,n_z,N_j,d_grid_cc,a_grid,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initialGE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_finalGE,n_d_cc,n_a,n_z,N_j,pi_z,ParamsGE,simoptions);

    PricePathGE.r=p_eqm.r*ones(1,T);
    PricePathGE.w=p_eqm.w*ones(1,T);
    ParamPathGE.sigma=Params.sigma*ones(1,T);
    Tbump=floor(2*T/3);
    PricePathBumped=PricePathGE;
    PricePathBumped.r(1:Tbump)=PricePathGE.r(1:Tbump).*(1+0.01+0.04*linspace(0,1,Tbump));
    PricePathBumped.w(1:Tbump)=PricePathGE.w(1:Tbump).*(1+0.05-0.04*linspace(0,1,Tbump));

    transpathoptionsS=transpathoptionsbaseline;
    transpathoptionsS.GEnewprice=3;
    transpathoptionsS.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
    transpathoptionsS.maxiter=250;
    transpathoptionsS.verbose=0;
    % Set here, not left to the default, because the checks in (3) read them. The distance below is
    % the one the solver uses under multiGEcriterion=1: the root sum of squares over the conditions in
    % each period, then the max over periods.
    transpathoptionsS.toleranceGEcondns=1e-4;
    transpathoptionsS.multiGEcriterion=1;

    %% (1) plain shooting
    tic;
    [PricePathB,GEcondnPathB]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d_cc, n_a, n_z, N_j, d_grid_cc,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptions);
    timeB=toc;

    %% (2) shooting with local search
    vfoptionsL=vfoptions;
    vfoptionsL.localsearch=1;
    vfoptionsL.nlocalsearch=1;
    tic;
    [PricePathL,GEcondnPathL]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d_cc, n_a, n_z, N_j, d_grid_cc,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptionsL);
    timeL=toc;
    fprintf('%s: plain shooting took %.2f seconds, with local search %.2f seconds \n',lab,timeB,timeL)
    % Both are solved to toleranceGEcondns from the same start, so they agree to about that tolerance
    % divided by the update factor, not exactly
    fprintf('%s: (2) price path with local search against plain shooting, r, should be close to zero: %.3e \n',lab,max(abs(PricePathL.r-PricePathB.r)))
    fprintf('%s: (2) price path with local search against plain shooting, w, should be close to zero: %.3e \n',lab,max(abs(PricePathL.w-PricePathB.w)))

    %% (3) CERTIFIED: the conditions returned come from the final, standard, sweep
    distB=max(sqrt(GEcondnPathB.CapitalMarket.^2+GEcondnPathB.LabourMarket.^2));
    distL=max(sqrt(GEcondnPathL.CapitalMarket.^2+GEcondnPathL.LabourMarket.^2));
    fprintf('%s: (3) plain shooting general eqm conditions within tolerance, should be one: %i \n',lab,distB<=transpathoptionsS.toleranceGEcondns)
    fprintf('%s: (3) local search general eqm conditions within tolerance, should be one: %i \n',lab,distL<=transpathoptionsS.toleranceGEcondns)
    % And independently of the solver's bookkeeping: one plain shooting iteration AT the local search
    % answer, which is converged at once if that answer really is the path
    transpathoptionsX=transpathoptionsS;
    transpathoptionsX.maxiter=1;
    [~,GEcondnPathX]=TransitionPath_Case1_FHorz(PricePathL, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d_cc, n_a, n_z, N_j, d_grid_cc,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsX, simoptions, vfoptions);
    distX=max(sqrt(GEcondnPathX.CapitalMarket.^2+GEcondnPathX.LabourMarket.^2));
    fprintf('%s: (3) the local search answer re-evaluated by an independent plain sweep, general eqm conditions within tolerance, should be one: %i \n',lab,distX<=transpathoptionsS.toleranceGEcondns)
    fprintf('%s: (3) local search conditions against the independent plain sweep at the same prices, this should be zero: %.3e \n',lab,max(abs(GEcondnPathL.CapitalMarket-GEcondnPathX.CapitalMarket))+max(abs(GEcondnPathL.LabourMarket-GEcondnPathX.LabourMarket)))
end

%% (4) and (5) once, on the with-d model without GI (the last set-up above is with-d WITH GI, so rebuild)
ReturnFn=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
FnsToEvaluateGE=struct();
FnsToEvaluateGE.K=@(d,aprime,a,z) a;
FnsToEvaluateGE.N=@(d,aprime,a,z,kappa_j) kappa_j*z*d;
vfoptions=struct();
simoptions=struct();
[V_finalGE,Policy_finalGE]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions);
AgentDist_initialGE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_finalGE,n_d,n_a,n_z,N_j,pi_z,ParamsGE,simoptions);
vfoptionsL=struct();
vfoptionsL.localsearch=1;
vfoptionsL.nlocalsearch=1;
lab='localsearch FHorz shooting';

%% (4) ValueFnOnTransPath_Case1_FHorz ignores localsearch=1, with a warning, and refuses 2
lastwarn('');
[VPath0,PolicyPath0]=ValueFnOnTransPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, Policy_finalGE, ParamsGE, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
[VPath1,PolicyPath1]=ValueFnOnTransPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, Policy_finalGE, ParamsGE, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsL);
[warnmsg,~]=lastwarn;
fprintf('%s: (4) ValueFnOnTransPath with localsearch=1 warns that it is ignored, should be one: %i \n',lab,contains(warnmsg,'localsearch=1 is IGNORED'))
fprintf('%s: (4) ValueFnOnTransPath with localsearch=1 against without, V, this should be zero: %.3e \n',lab,max(abs(VPath1(:)-VPath0(:))))
fprintf('%s: (4) ValueFnOnTransPath with localsearch=1 against without, Policy, this should be zero: %.3e \n',lab,max(abs(PolicyPath1(:)-PolicyPath0(:))))
clear VPath0 VPath1 PolicyPath0 PolicyPath1
err4=0;
vfoptionsE=vfoptionsL; vfoptionsE.localsearch=2;
try
    ValueFnOnTransPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, Policy_finalGE, ParamsGE, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsE);
catch ME
    err4=contains(ME.message,'not yet implemented for finite horizon');
end
fprintf('%s: (4) ValueFnOnTransPath with localsearch=2 is refused by name, should be one: %i \n',lab,err4)

%% (5) The refusals in TransitionPath_Case1_FHorz, as a set
% Every one of these runs BEFORE any solving, so each costs nothing. maxiter=1 all the same, so that a
% refusal that has been lifted costs one iteration rather than a full solve before it shows up as a 0.
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
    fired=contains(ME.message,'requires transpathoptions.GEnewprice=3');
end
refusalname{end+1}='GEnewprice=1'; refusalfired(end+1)=fired;
% Anderson, not yet wired up for FHorz
transpathoptionsE=rmfield(transpathoptionsS,'GEnewprice3'); transpathoptionsE.GEnewprice=2; transpathoptionsE.GEnewprice2.howtoupdate=transpathoptionsS.GEnewprice3.howtoupdate; fired=0;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptionsL);
catch ME
    fired=contains(ME.message,'requires transpathoptions.GEnewprice=3');
end
refusalname{end+1}='GEnewprice=2'; refusalfired(end+1)=fired;
% fastOLG=1
transpathoptionsE=transpathoptionsS; transpathoptionsE.fastOLG=1; fired=0;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptionsL);
catch ME
    fired=contains(ME.message,'not yet implemented with transpathoptions.fastOLG=1');
end
refusalname{end+1}='fastOLG=1'; refusalfired(end+1)=fired;
% divide-and-conquer
vfoptionsE=vfoptionsL; vfoptionsE.divideandconquer=1; fired=0;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptionsE);
catch ME
    fired=contains(ME.message,'together with vfoptions.divideandconquer=1');
end
refusalname{end+1}='divideandconquer=1'; refusalfired(end+1)=fired;
% an e variable
vfoptionsE=vfoptionsL; vfoptionsE.n_e=3; vfoptionsE.e_grid=linspace(0.9,1.1,3)'; vfoptionsE.pi_e=ones(3,1)/3; fired=0;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptionsE);
catch ME
    fired=contains(ME.message,'not implemented with an e variable');
end
refusalname{end+1}='e variable'; refusalfired(end+1)=fired;
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
