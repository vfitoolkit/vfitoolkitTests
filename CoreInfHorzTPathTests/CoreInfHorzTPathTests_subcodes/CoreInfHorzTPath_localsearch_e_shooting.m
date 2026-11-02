function output=CoreInfHorzTPath_localsearch_e_shooting(T,n_d_GE,n_a_GE,n_z,d_grid_GE,a_grid_GE,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline)
% vfoptions.localsearch=1 with z AND an iid e, inside a general eqm transition path: section (8) of stage L1 of
% InfHorzTPath_localsearch_noz_e_proposal.md. The value fn step and seed-sweep sections, (1)-(7) and (9), are
% CoreInfHorzTPath_localsearch_e (doPart(6)); this is the general eqm half, split out of it on 2026-11-01 when doPart(11)
% was merged into doPart(6) and doPart(7). Read that header for the reference per (d,a,z,e).

gpudev=gpuDevice;

n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;

%% (8) Inside a general eqm transition path, without and with d
% The model of doPart(7) with an e: K is assets, N is efficiency labour (z*e, or d*z*e with d). Shooting and Anderson,
% each standard and with local search; then the local search answer is re-evaluated with one exact sweep (zero
% damping, maxiter=1), which is the certificate.
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;
for dcase=1:2
    vfoptions=struct();
    vfoptions.n_e=n_e;
    vfoptions.e_grid=e_grid;
    vfoptions.pi_e=pi_e;
    simoptions=struct();
    simoptions.n_e=n_e;
    simoptions.e_grid=e_grid;
    simoptions.pi_e=pi_e;
    FnsToEvaluateGE=struct();
    if dcase==1
        lab='localsearch e GE';
        n_dc=0; d_gridc=[];
        ReturnFn=@(aprime,a,z,e,r,w,sigma) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,sigma);
        FnsToEvaluateGE.K=@(aprime,a,z,e) a;
        FnsToEvaluateGE.N=@(aprime,a,z,e) z*e;
    else
        lab='localsearch d e GE';
        n_dc=n_d_GE; d_gridc=d_grid_GE;
        ReturnFn=@(d,aprime,a,z,e,r,w,sigma,eta,varphi) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,sigma,eta,varphi);
        FnsToEvaluateGE.K=@(d,aprime,a,z,e) a;
        FnsToEvaluateGE.N=@(d,aprime,a,z,e) d*z*e;
    end
    GeneralEqmEqnsGE=struct();
    GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
    GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

    Params.r=0.04; Params.w=1; % initial guess
    [p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_dc, n_a_GE, n_z, 0, pi_z, d_gridc, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, struct(), simoptions, vfoptions);
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
