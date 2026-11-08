function output=EZFHorzTPath_fastOLG2(T,n_d,n_a,n_z,n_e,N_j,d_grid,a_grid,z_grid,pi_z,e_grid,pi_e,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline)
% transpathoptions.fastOLG=2 with Epstein-Zin preferences: fastOLG for the local search sweeps, slowOLG for the
% full sweeps (the first of a solve, and every seed or verification sweep). No figures.
% The Epstein-Zin mirror of CoreFHorzTPath_fastOLG2 (CoreFHorzTPathTests doPart(14)). The model is the one in
% EZFHorzTPath_localsearch_shooting (doPart(7)): a Cobb-Douglas firm supplies both r and w, and the path is a
% null reform solved from a bumped starting guess. Survival probabilities and the warm-glow are on throughout.
% The rows are the fastOLG rows of doPart(7), which converge there.
%
% One solve per row of the configuration table below, done twice with local search, by the SAME algorithm
% with the SAME options: once with fastOLG=1, once with fastOLG=2.
%   (0) the stationary general eqm the path ends at is within the path's toleranceGEcondns
%   (1) the two give the same price path and the same general eqm conditions, exactly: the full sweeps
%       differ only in layout, and the slowOLG and fastOLG Epstein-Zin raws solve the same problem to the
%       same answer (doParts 1-2 check that, sweep by sweep)
%   (2) the DISCRIMINATOR for (1): the fastOLG=2 solve did solve its full sweeps age by age. Run with
%       transpathoptions.verbose=1, each such sweep prints one line; the count must be at least one
%       (iteration 1 is always a full sweep), or (1) compared fastOLG=1 with itself
% Then, once, on the first row's model and price path:
%   (3) ValueFnOnTransPath_Case1_FHorz: its single pass is a full sweep, so fastOLG=2 must give exactly
%       fastOLG=0; with local search forced (vfoptions.localsearch=2) exactly fastOLG=1

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

% columns: EZ case (1=cons units, 2=positive utils, 3=negative utils), shocks (1=z,2=noz,3=e,4=noz_e), hasd, GI, GEnewprice, divideandconquer
config=[1,1,0,0,3,0;   % shooting, the three EZ cases
        3,1,0,1,3,0;
        2,1,1,0,3,0;   % with d
        1,1,0,1,3,1;   % divide-and-conquer (honoured on the full sweeps)
        2,3,0,1,3,0;   % z and e
        3,2,1,1,3,0;   % no z
        3,4,0,0,3,0;   % e without z
        2,1,0,1,2,0];  % Anderson
casename={'cons units','positive utils','negative utils'};
casefn={'cons','positiveUtils','negativeUtils'};
shockname={'z','noz','e','noz_e'};
algoname={'','Anderson','shooting'};

for cc=1:size(config,1)
    ezcase=config(cc,1); ss=config(cc,2); hasd=config(cc,3); useGI=config(cc,4); GEnewprice=config(cc,5); useDC=config(cc,6);
    hasz=(ss==1 || ss==3);
    hase=(ss==3 || ss==4);
    lab=sprintf('EZ fastOLG=2 %s (%s, %s, d%i, GI%i, DC%i)',algoname{GEnewprice},casename{ezcase},shockname{ss},hasd,useGI,useDC);

    vfoptions=struct();
    simoptions=struct();
    vfoptions.exoticpreferences='EpsteinZin';
    if ezcase==1
        vfoptions.EZutils=0; vfoptions.EZriskaversion='ezgamma'; vfoptions.EZeis='ezphi';
        vfoptions.WarmGlowBequestsFn=@(aprime,wg1,wg2) EZWarmGlowFn_cons(aprime,wg1,wg2);
    elseif ezcase==2
        vfoptions.EZutils=1; vfoptions.EZpositiveutility=1; vfoptions.EZriskaversion='ezrisk';
        vfoptions.WarmGlowBequestsFn=@(aprime,wg1,wg2,wg3) EZWarmGlowFn_positiveUtils(aprime,wg1,wg2,wg3);
    else
        vfoptions.EZutils=1; vfoptions.EZpositiveutility=0; vfoptions.EZriskaversion='ezrisk';
        vfoptions.WarmGlowBequestsFn=@(aprime,wg1,wg2,wg3) EZWarmGlowFn_negativeUtils(aprime,wg1,wg2,wg3);
    end
    vfoptions.survivalprobability='sj';
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
    % Return fn for this model, from the bank's EZ set: the states, then the parameters
    statenames={};
    if hasd
        statenames{end+1}='d';
    end
    statenames=[statenames,{'aprime','a'}];
    if hasz
        statenames{end+1}='z';
    end
    if hase
        statenames{end+1}='e';
    end
    ReturnFnParamNames={'r','w','kappa_j'};
    if ezcase>1
        ReturnFnParamNames{end+1}='ezsigma';
    end
    if hasd
        ReturnFnParamNames{end+1}='varphi';
    end
    ReturnFnParamNames=[ReturnFnParamNames,{'agej','Jr','pension'}];
    if hasd
        dname='d';
    else
        dname='nod';
    end
    if hasz
        zname='z';
    else
        zname='noz';
    end
    if hase
        ename='e';
    else
        ename='noe';
    end
    arglist=strjoin([statenames,ReturnFnParamNames],',');
    ReturnFn=str2func(['@(',arglist,') EZReturnFn_',casefn{ezcase},'_',dname,'_',zname,'_',ename,'_nosemiz(',arglist,')']);
    % Aggregates: K is assets, N is efficiency labour, kappa_j times the shocks (times d)
    stateargs=strjoin(statenames,',');
    FnsToEvaluateGE=struct();
    FnsToEvaluateGE.K=str2func(['@(',stateargs,') a']);
    Nexpr='kappa_j';
    if hasz
        Nexpr=[Nexpr,'*z'];
    end
    if hase
        Nexpr=[Nexpr,'*e'];
    end
    if hasd
        Nexpr=[Nexpr,'*d'];
    end
    FnsToEvaluateGE.N=str2func(['@(',stateargs,',kappa_j) ',Nexpr]);
    GeneralEqmEqnsGE=struct();
    GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
    GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

    %% Stationary general eqm, then a null-reform path from a bumped starting guess
    % As in doPart(7): solved to 1e-8, so that it sits inside the path's toleranceGEcondns (1e-4)
    heteroagentoptionsGE=struct();
    heteroagentoptionsGE.verbose=0;
    heteroagentoptionsGE.toleranceGEprices=1e-8;
    heteroagentoptionsGE.toleranceGEcondns=1e-8;
    heteroagentoptionsGE.maxiter=2000;
    heteroagentoptionsGE.outputGEstruct=2;
    [p_eqm,GEcondnsStat]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d_cc,n_a,n_z_cc,N_j,0,pi_z_cc,d_grid_cc,a_grid,z_grid_cc,ReturnFn,FnsToEvaluateGE,GeneralEqmEqnsGE,Params,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptionsGE,simoptions,vfoptions);
    % (0) root sum of squares, the same distance the path uses per period under multiGEcriterion=1
    distStat=sqrt(sum(gather(GEcondnsStat(:)).^2));
    fprintf('%s: (0) stationary general eqm conditions, root sum of squares, within the path tolerance, should be one: %i (it is %8.3e) \n',lab,distStat<=1e-4,distStat) % 1e-4 is the path's toleranceGEcondns, set below
    ParamsGE=Params;
    ParamsGE.r=p_eqm.r;
    ParamsGE.w=p_eqm.w;
    [V_finalGE,Policy_finalGE]=ValueFnIter_Case1_FHorz(n_d_cc,n_a,n_z_cc,N_j,d_grid_cc,a_grid,z_grid_cc,pi_z_cc,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initialGE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_finalGE,n_d_cc,n_a,n_z_cc,N_j,pi_z_cc,ParamsGE,simoptions);

    PricePathGE.r=p_eqm.r*ones(1,T);
    PricePathGE.w=p_eqm.w*ones(1,T);
    ParamPathGE=struct();
    ParamPathGE.firmA=Params.firmA*ones(1,T);
    Tbump=floor(2*T/3);
    PricePathBumped=PricePathGE;
    PricePathBumped.r(1:Tbump)=PricePathGE.r(1:Tbump).*(1+0.01+0.04*linspace(0,1,Tbump));
    PricePathBumped.w(1:Tbump)=PricePathGE.w(1:Tbump).*(1+0.05-0.04*linspace(0,1,Tbump));

    transpathoptionsS=transpathoptionsbaseline;
    transpathoptionsS.GEnewprice=GEnewprice;
    if GEnewprice==3
        transpathoptionsS.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
    else
        transpathoptionsS.GEnewprice2.howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
    end
    transpathoptionsS.maxiter=1000;
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
    tic;
    [PricePath1,GEcondnPath1]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d_cc, n_a, n_z_cc, N_j, d_grid_cc,a_grid,z_grid_cc, pi_z_cc, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions1, simoptions, vfoptionsL);
    time1=toc;
    transpathoptions2=transpathoptionsS;
    transpathoptions2.fastOLG=2;
    transpathoptions2.verbose=1; % so each full sweep solved age by age says so, see (2)
    tic;
    out2=evalc('[PricePath2,GEcondnPath2]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d_cc, n_a, n_z_cc, N_j, d_grid_cc,a_grid,z_grid_cc, pi_z_cc, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions2, simoptions, vfoptionsL);');
    time2=toc;
    fprintf('%s: fastOLG=1 took %.2f seconds, fastOLG=2 %.2f seconds \n',lab,time1,time2)
    fprintf('%s: (1) price path against fastOLG=1, r, this should be zero: %.3e \n',lab,max(abs(PricePath2.r-PricePath1.r)))
    fprintf('%s: (1) price path against fastOLG=1, w, this should be zero: %.3e \n',lab,max(abs(PricePath2.w-PricePath1.w)))
    fprintf('%s: (1) general eqm conditions against fastOLG=1, this should be zero: %.3e \n',lab,max(abs(GEcondnPath2.CapitalMarket-GEcondnPath1.CapitalMarket))+max(abs(GEcondnPath2.LabourMarket-GEcondnPath1.LabourMarket)))
    dist1=max(sqrt(GEcondnPath1.CapitalMarket.^2+GEcondnPath1.LabourMarket.^2));
    fprintf('%s: (1) fastOLG=1 general eqm conditions within tolerance, should be one: %i \n',lab,dist1<=transpathoptionsS.toleranceGEcondns)
    % anti-vacuity: the solve travelled from the bumped guess
    fprintf('%s: (1) the solve moved the price path off its starting guess, should NOT be zero: %.3e \n',lab,max(abs(PricePath1.r-PricePathBumped.r)))
    %% (2) the discriminator
    nfull=numel(strfind(out2,'a full sweep, solved age by age'));
    fprintf('%s: (2) full sweeps solved age by age (one line each, under verbose), should NOT be zero: %i \n',lab,nfull)

    %% (3) once, on the first row: ValueFnOnTransPath_Case1_FHorz with fastOLG=2
    if cc==1
        lab3='EZ fastOLG=2 ValueFnOnTransPath (cons units, z)';
        VPathL=cell(1,3); PolicyPathL=cell(1,3);
        for fo=0:2
            transpathoptionsV=transpathoptionsbaseline; transpathoptionsV.fastOLG=fo;
            [VPathL{fo+1},PolicyPathL{fo+1}]=ValueFnOnTransPath_Case1_FHorz(PricePath1, ParamPathGE, T, V_finalGE, Policy_finalGE, ParamsGE, n_d_cc, n_a, n_z_cc, N_j, d_grid_cc, a_grid,z_grid_cc, pi_z_cc, DiscountFactorParamNames, ReturnFn, transpathoptionsV, vfoptions);
        end
        fprintf('%s: (3) the single pass is a full sweep, against fastOLG=0, this should be zero, V: %.3e \n',lab3,max(abs(VPathL{3}(:)-VPathL{1}(:))))
        fprintf('%s: (3) the single pass is a full sweep, against fastOLG=0, this should be zero, Policy: %.3e \n',lab3,max(abs(PolicyPathL{3}(:)-PolicyPathL{1}(:))))
        vfoptionsF=vfoptions; vfoptionsF.localsearch=2; vfoptionsF.nlocalsearch=3;
        for fo=1:2
            transpathoptionsV=transpathoptionsbaseline; transpathoptionsV.fastOLG=fo;
            [VPathL{fo+1},PolicyPathL{fo+1}]=ValueFnOnTransPath_Case1_FHorz(PricePath1, ParamPathGE, T, V_finalGE, Policy_finalGE, ParamsGE, n_d_cc, n_a, n_z_cc, N_j, d_grid_cc, a_grid,z_grid_cc, pi_z_cc, DiscountFactorParamNames, ReturnFn, transpathoptionsV, vfoptionsF);
        end
        fprintf('%s: (3) with local search forced, against fastOLG=1, this should be zero, V: %.3e \n',lab3,max(abs(VPathL{3}(:)-VPathL{2}(:))))
        fprintf('%s: (3) with local search forced, against fastOLG=1, this should be zero, Policy: %.3e \n',lab3,max(abs(PolicyPathL{3}(:)-PolicyPathL{2}(:))))
        clear VPathL PolicyPathL
    end
end

output=1;

end
