function output=CoreInfHorzTPath_localsearch_noz_shooting(T,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,Params,DiscountFactorParamNames,transpathoptionsbaseline)
% vfoptions.localsearch=1 without z, with an iid e, inside a general eqm transition path: section (8) of stages L2/L3 of
% InfHorzTPath_localsearch_noz_e_proposal.md. The value fn step and seed-sweep sections are CoreInfHorzTPath_localsearch_noz
% (doPart(6)); this is the general eqm half, split out of it on 2026-11-01 when doPart(11) was merged into doPart(6) and
% doPart(7). Without z and without e the stationary general eqm is degenerate, so this is run with e only.

gpudev=gpuDevice;

n_z=0; z_grid=[]; pi_z=[]; % no z: the inputs stay positional, as empty

%% (8) Inside a general eqm transition path, no z with e, without and with d
% The model of CoreInfHorzTPath_localsearch_e section (8) without z: K is assets, N is efficiency labour (e, or d*e).
% The stationary solve is started from a GI solve, as in CoreInfHorzTPath_d_noz_e_nosemiz (K and N are step functions
% of the prices on coarse grids and the default start stalls).
% THE e HERE IS RISKIER THAN THE BANK'S (5 points, sd 0.4, against 3 points, sd 0.1). Without z, the bank's e is
% almost no risk: there is next to no precautionary saving, r sits at 1/beta-1 where asset demand is nearly flat, and
% the stationary solve stalled (CapitalMarket residual 1.7e-2 without d and 5.4e-2 with d on 2026-11-01), which made
% toleranceGEcondns=5*residual so loose that every check below passed without testing anything (Anderson with and
% without local search agreed to 0, as both stopped at the first evaluation). Hence also the guard: the stationary
% eqm must be solved to 1e-3, or the section says so as a failure rather than passing vacuously.
[e_gridGE,pi_eGE]=discretizeAR1_FarmerToda(0,0,0.4,5);
pi_eGE=pi_eGE(1,:)';
e_gridGE=exp(e_gridGE);
n_eGE=5;
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;
for dcase=1:2
    vfoptions=struct();
    vfoptions.n_e=n_eGE;
    vfoptions.e_grid=e_gridGE;
    vfoptions.pi_e=pi_eGE;
    simoptions=struct();
    simoptions.n_e=n_eGE;
    simoptions.e_grid=e_gridGE;
    simoptions.pi_e=pi_eGE;
    FnsToEvaluateGE=struct();
    if dcase==1
        lab='localsearch noz e GE';
        n_dc=0; d_gridc=[];
        ReturnFn=@(aprime,a,e,r,w,sigma) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,sigma);
        FnsToEvaluateGE.K=@(aprime,a,e) a;
        FnsToEvaluateGE.N=@(aprime,a,e) e;
    else
        lab='localsearch d noz e GE';
        n_dc=n_d_GE; d_gridc=d_grid_GE;
        ReturnFn=@(d,aprime,a,e,r,w,sigma,eta,varphi) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,sigma,eta,varphi);
        FnsToEvaluateGE.K=@(d,aprime,a,e) a;
        FnsToEvaluateGE.N=@(d,aprime,a,e) d*e;
    end
    GeneralEqmEqnsGE=struct();
    GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
    GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

    Params.r=0.04; Params.w=1; % initial guess
    ParamsStart=Params;
    % Both models get a start from a GI solve: without d as well, on this grid the default start stalled
    % (CapitalMarket residual 3.2e-3 on 2026-11-01, caught by the guard below); with d it got 3.1e-4 from the GI start
    vfoptionsGIs=vfoptions; vfoptionsGIs.gridinterplayer=1; vfoptionsGIs.ngridinterp=5;
    simoptionsGIs=simoptions; simoptionsGIs.gridinterplayer=1; simoptionsGIs.ngridinterp=5;
    p_eqmstart=HeteroAgentStationaryEqm_InfHorz(n_dc, n_a_GE, n_z, 0, pi_z, d_gridc, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, struct(), simoptionsGIs, vfoptionsGIs);
    ParamsStart.r=p_eqmstart.r;
    ParamsStart.w=p_eqmstart.w;
    [p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_dc, n_a_GE, n_z, 0, pi_z, d_gridc, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsStart, DiscountFactorParamNames, [], [], [], {'r','w'}, struct(), simoptions, vfoptions);
    fprintf('%s: stationary GE r=%2.8f, w=%2.8f, conditions CapitalMarket=%.3e LabourMarket=%.3e \n',lab,p_eqm.r,p_eqm.w,GEcondns.CapitalMarket,GEcondns.LabourMarket)
    ParamsGE=Params;
    ParamsGE.r=p_eqm.r;
    ParamsGE.w=p_eqm.w;
    [V_finalGE,Policy_finalGE]=ValueFnIter_InfHorz(n_dc,n_a_GE,n_z,d_gridc,a_grid_GE,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initialGE=StationaryDist_InfHorz(Policy_finalGE,n_dc,n_a_GE,n_z,pi_z,simoptions,ParamsGE,[]);

    PricePathGE.r=p_eqm.r*ones(1,T);
    PricePathGE.w=p_eqm.w*ones(1,T);
    ParamPathGE=struct();
    ParamPathGE.sigma=Params.sigma*ones(1,T); % constant: a null reform, so the answer is the stationary eqm
    Tbump=floor(2*T/3);
    PricePathBumped=PricePathGE;
    PricePathBumped.r(1:Tbump)=PricePathGE.r(1:Tbump).*(1+0.01+0.04*linspace(0,1,Tbump));
    PricePathBumped.w(1:Tbump)=PricePathGE.w(1:Tbump).*(1+0.05-0.04*linspace(0,1,Tbump));

    GEfloor=max(abs(GEcondns.CapitalMarket),abs(GEcondns.LabourMarket));
    fprintf('%s: (8) the stationary eqm is solved to 1e-3, so the checks below are not vacuous, should be one: %i \n',lab,GEfloor<=1e-3)
    transpathoptionsS=transpathoptionsbaseline;
    transpathoptionsS.GEnewprice=3;
    transpathoptionsS.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
    transpathoptionsS.maxiter=1000;
    transpathoptionsS.verbose=0;
    transpathoptionsS.toleranceGEprices=10^(-4);
    transpathoptionsS.toleranceGEcondns=5*GEfloor;
    transpathoptionsAnd=transpathoptionsbaseline;
    transpathoptionsAnd.GEnewprice=2;
    transpathoptionsAnd.GEnewprice2.howtoupdate=transpathoptionsS.GEnewprice3.howtoupdate;
    transpathoptionsAnd.maxiter=200;
    transpathoptionsAnd.verbose=0;
    transpathoptionsAnd.toleranceGEprices=transpathoptionsS.toleranceGEprices;
    transpathoptionsAnd.toleranceGEcondns=transpathoptionsS.toleranceGEcondns;

    vfoptionsLS=vfoptions;
    vfoptionsLS.localsearch=1;
    vfoptionsLS.nlocalsearch=1;

    tic;
    PricePathA=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_dc, n_a_GE, n_z, d_gridc,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsS, simoptions, vfoptions, []);
    wait(gpudev); tA=toc;
    tic;
    PricePathB=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_dc, n_a_GE, n_z, d_gridc,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsS, simoptions, vfoptionsLS, []);
    wait(gpudev); tB=toc;
    PricePathAnd=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_dc, n_a_GE, n_z, d_gridc,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsAnd, simoptions, vfoptions, []);
    PricePathAndLS=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_dc, n_a_GE, n_z, d_gridc,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsAnd, simoptions, vfoptionsLS, []);
    fprintf('%s: for information, standard shooting %6.2fs, shooting with local search from nlocalsearch=1 %6.2fs \n',lab,tA,tB)
    fprintf('%s: (8) shooting with local search vs the standard shooting solve, r, this should be close to zero: %.3e \n',lab,max(abs(gather(PricePathB.r)-gather(PricePathA.r))))
    fprintf('%s: (8) same, w, this should be close to zero: %.3e \n',lab,max(abs(gather(PricePathB.w)-gather(PricePathA.w))))
    fprintf('%s: (8) Anderson with local search vs the standard Anderson solve, r, this should be close to zero: %.3e \n',lab,max(abs(gather(PricePathAndLS.r)-gather(PricePathAnd.r))))
    fprintf('%s: (8) same, w, this should be close to zero: %.3e \n',lab,max(abs(gather(PricePathAndLS.w)-gather(PricePathAnd.w))))
    fprintf('%s: (8) for information, the standard shooting solve deviates from the stationary eqm by at most r %.3e, w %.3e \n',lab,max(abs(gather(PricePathA.r)-p_eqm.r)),max(abs(gather(PricePathA.w)-p_eqm.w)))
    % The certificate: re-evaluate each local search answer with one exact sweep and zero damping
    transpathoptionsV=transpathoptionsbaseline;
    transpathoptionsV.GEnewprice=3;
    transpathoptionsV.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0; 'LabourMarket','w',0,0};
    transpathoptionsV.maxiter=1;
    transpathoptionsV.verbose=0;
    [~,GEcondnPathV]=TransitionPath_InfHorz(PricePathB, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_dc, n_a_GE, n_z, d_gridc,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsV, simoptions, vfoptions, []);
    GEdistV=max(max(abs(gather(GEcondnPathV.CapitalMarket))),max(abs(gather(GEcondnPathV.LabourMarket))));
    [~,GEcondnPathVA]=TransitionPath_InfHorz(PricePathAndLS, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_dc, n_a_GE, n_z, d_gridc,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsV, simoptions, vfoptions, []);
    GEdistVA=max(max(abs(gather(GEcondnPathVA.CapitalMarket))),max(abs(gather(GEcondnPathVA.LabourMarket))));
    fprintf('%s: (8) the shooting local search answer re-evaluated with an exact sweep is within toleranceGEcondns, should be one: %i \n',lab,GEdistV<=transpathoptionsS.toleranceGEcondns)
    fprintf('%s: (8) the Anderson local search answer re-evaluated with an exact sweep is within toleranceGEcondns, should be one: %i \n',lab,GEdistVA<=transpathoptionsS.toleranceGEcondns)
    fprintf('%s: (8) for information, those distances are %.3e and %.3e against toleranceGEcondns %.3e \n',lab,GEdistV,GEdistVA,transpathoptionsS.toleranceGEcondns)
end

output=struct();

end
