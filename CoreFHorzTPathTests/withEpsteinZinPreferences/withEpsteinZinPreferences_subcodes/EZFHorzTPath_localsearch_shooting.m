function output=EZFHorzTPath_localsearch_shooting(T,n_d,n_a,n_z,n_e,N_j,d_grid,a_grid,z_grid,pi_z,e_grid,pi_e,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline)
% vfoptions.localsearch=1 with Epstein-Zin preferences inside a real general eqm transition path. No figures.
% The Epstein-Zin mirror of CoreFHorzTPath_localsearch_shooting (CoreFHorzTPathTests doPart(10)): a Cobb-Douglas
% firm supplies both r and w, and the path is a null reform solved from a bumped starting guess, so the solver
% has to travel. Survival probabilities and the warm-glow of bequests are on throughout.
%
% One solve per row of the configuration table below: each is solved twice, once plainly and once with
% local search, by the SAME algorithm on the SAME layout with the SAME divideandconquer setting:
%   (2) the two must converge to the same price path
%   (3) the local search answer is CERTIFIED: its returned general eqm conditions are within tolerance,
%       and an INDEPENDENT plain sweep at its prices gives the same conditions
% The table covers the three Epstein-Zin cases, both layouts, shooting and Anderson, divide-and-conquer
% (honoured on the seed sweeps), and the four shock cases, each with or without d and GI.
% Then, once, (4): ValueFnOnTransPath_Case1_FHorz with Epstein-Zin ignores localsearch=1 with a warning, and
% localsearch=2 forces it (at the cap it must equal the plain value fn).

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

% columns: EZ case (1=cons units, 2=positive utils, 3=negative utils), fastOLG, shocks (1=z,2=noz,3=e,4=noz_e), hasd, GI, GEnewprice, divideandconquer
config=[1,0,1,0,0,3,0;   % the three EZ cases, slowOLG, z, by shooting
        2,0,1,0,0,3,0;   % rows 2 and 14 use n_a=301, rows 4, 8 and 10 n_d=21 (see below)
        3,0,1,0,0,3,0;
        3,0,1,1,1,3,0;   % with d and GI
        1,1,1,0,0,3,0;   % fastOLG
        3,1,1,0,1,3,0;
        2,1,1,1,0,3,0;
        3,0,1,1,0,2,0;   % Anderson
        2,1,1,0,1,2,0;
        3,0,1,1,0,3,1;   % divide-and-conquer, honoured on the seed sweeps
        1,1,1,0,1,3,1;
        3,0,3,1,0,3,0;   % z and e
        2,1,3,0,1,3,0;
        1,0,2,0,0,3,0;   % no z (the certainty equivalent is the identity, but the composition is not)
        3,1,2,1,1,3,0;
        2,0,4,1,1,3,0;   % e without z
        3,1,4,0,0,3,0];
casename={'cons units','positive utils','negative utils'};
casefn={'cons','positiveUtils','negativeUtils'};
layoutname={'slowOLG','fastOLG'};
shockname={'z','noz','e','noz_e'};
algoname={'','Anderson','shooting'};

for cc=1:size(config,1)
    ezcase=config(cc,1); fastOLG=config(cc,2); ss=config(cc,3); hasd=config(cc,4); useGI=config(cc,5); GEnewprice=config(cc,6); useDC=config(cc,7);
    hasz=(ss==1 || ss==3);
    hase=(ss==3 || ss==4);
    lab=sprintf('EZ localsearch FHorz %s (%s, %s, %s, d%i, GI%i, DC%i)',algoname{GEnewprice},casename{ezcase},layoutname{fastOLG+1},shockname{ss},hasd,useGI,useDC);

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
        if cc==4 || cc==8 || cc==10
            % negative utils, slowOLG, z, with d: at n_d=9 (steps of 0.125) neither the stationary general eqm nor
            % the path got within tolerance, with or without GI, so a finer d grid on the same [0,1]
            n_d_cc=21; d_grid_cc=linspace(0,1,n_d_cc)';
        end
    else
        n_d_cc=0; d_grid_cc=[];
    end
    n_a_cc=n_a; a_grid_cc=a_grid;
    if cc==2 || cc==14
        % without GI at n_a=101 K jumps between grid steps (positive utils: 0.2507/0.2531 at r=0.2533, w=0.2384),
        % so no price gets the conditions within 1e-4; a finer grid of the same shape keeps these rows GI0
        n_a_cc=301; a_grid_cc=5*linspace(0,1,n_a_cc)'.^3;
    end
    % zero assets, mid points for the shocks
    if ss==1
        jequaloneDist=zeros([n_a_cc,n_z],'gpuArray'); jequaloneDist(1,ceil(n_z/2))=1;
    elseif ss==2
        jequaloneDist=zeros([n_a_cc,1],'gpuArray'); jequaloneDist(1)=1;
    elseif ss==3
        jequaloneDist=zeros([n_a_cc,n_z,n_e],'gpuArray'); jequaloneDist(1,ceil(n_z/2),ceil(n_e/2))=1;
    else
        jequaloneDist=zeros([n_a_cc,n_e],'gpuArray'); jequaloneDist(1,ceil(n_e/2))=1;
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
    % The stationary solve and V_final use the same vfoptions WITHOUT local search (they are not path solves)
    % The path ends at this stationary general eqm, so it must sit well inside the path's toleranceGEcondns
    % (1e-4, set below), or the tail of the path cannot get within it and the path runs to maxiter. At the
    % default 1e-4 fminsearch stopped at 3.2e-4 for positive utils. outputGEstruct=2 returns the conditions
    % as a vector, for the check (0) below.
    heteroagentoptionsGE=struct();
    heteroagentoptionsGE.verbose=0;
    heteroagentoptionsGE.toleranceGEprices=1e-8;
    heteroagentoptionsGE.toleranceGEcondns=1e-8;
    heteroagentoptionsGE.maxiter=2000;
    heteroagentoptionsGE.outputGEstruct=2;
    [p_eqm,GEcondnsStat]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d_cc,n_a_cc,n_z_cc,N_j,0,pi_z_cc,d_grid_cc,a_grid_cc,z_grid_cc,ReturnFn,FnsToEvaluateGE,GeneralEqmEqnsGE,Params,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptionsGE,simoptions,vfoptions);
    % (0) root sum of squares, the same distance the path uses per period under multiGEcriterion=1
    distStat=sqrt(sum(gather(GEcondnsStat(:)).^2));
    fprintf('%s: (0) stationary general eqm conditions, root sum of squares, within the path tolerance, should be one: %i (it is %8.3e) \n',lab,distStat<=1e-4,distStat) % 1e-4 is the path's toleranceGEcondns, set below
    ParamsGE=Params;
    ParamsGE.r=p_eqm.r;
    ParamsGE.w=p_eqm.w;
    [V_finalGE,Policy_finalGE]=ValueFnIter_Case1_FHorz(n_d_cc,n_a_cc,n_z_cc,N_j,d_grid_cc,a_grid_cc,z_grid_cc,pi_z_cc,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initialGE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_finalGE,n_d_cc,n_a_cc,n_z_cc,N_j,pi_z_cc,ParamsGE,simoptions);

    PricePathGE.r=p_eqm.r*ones(1,T);
    PricePathGE.w=p_eqm.w*ones(1,T);
    ParamPathGE=struct();
    ParamPathGE.firmA=Params.firmA*ones(1,T);
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
    transpathoptionsS.maxiter=1000;
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
    [PricePathB,GEcondnPathB]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d_cc, n_a_cc, n_z_cc, N_j, d_grid_cc,a_grid_cc,z_grid_cc, pi_z_cc, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptions);
    timeB=toc;

    %% (2) the same solve with local search
    vfoptionsL=vfoptions;
    vfoptionsL.localsearch=1;
    vfoptionsL.nlocalsearch=1;
    tic;
    [PricePathL,GEcondnPathL]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d_cc, n_a_cc, n_z_cc, N_j, d_grid_cc,a_grid_cc,z_grid_cc, pi_z_cc, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptionsL);
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
    [~,GEcondnPathX]=TransitionPath_Case1_FHorz(PricePathL, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d_cc, n_a_cc, n_z_cc, N_j, d_grid_cc,a_grid_cc,z_grid_cc, pi_z_cc, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsX, simoptions, vfoptions);
    distX=max(sqrt(GEcondnPathX.CapitalMarket.^2+GEcondnPathX.LabourMarket.^2));
    fprintf('%s: (3) the local search answer re-evaluated by an independent plain sweep, general eqm conditions within tolerance, should be one: %i \n',lab,distX<=transpathoptionsS.toleranceGEcondns)
    if GEnewprice==3
        % Shooting returns the conditions of its final (standard) sweep at the returned prices, so these agree exactly
        fprintf('%s: (3) local search conditions against the independent plain sweep at the same prices, this should be zero: %.3e \n',lab,max(abs(GEcondnPathL.CapitalMarket-GEcondnPathX.CapitalMarket))+max(abs(GEcondnPathL.LabourMarket-GEcondnPathX.LabourMarket)))
    end

    %% (4) once, on the first row: ValueFnOnTransPath_Case1_FHorz with Epstein-Zin
    if cc==1
        lab4='EZ localsearch FHorz ValueFnOnTransPath';
        [VPath0,PolicyPath0]=ValueFnOnTransPath_Case1_FHorz(PricePathB, ParamPathGE, T, V_finalGE, Policy_finalGE, ParamsGE, n_d_cc, n_a_cc, n_z_cc, N_j, d_grid_cc, a_grid_cc,z_grid_cc, pi_z_cc, DiscountFactorParamNames, ReturnFn, transpathoptionsS, vfoptions);
        lastwarn('');
        [VPath1,PolicyPath1]=ValueFnOnTransPath_Case1_FHorz(PricePathB, ParamPathGE, T, V_finalGE, Policy_finalGE, ParamsGE, n_d_cc, n_a_cc, n_z_cc, N_j, d_grid_cc, a_grid_cc,z_grid_cc, pi_z_cc, DiscountFactorParamNames, ReturnFn, transpathoptionsS, vfoptionsL);
        warnmsg=lastwarn;
        fprintf('%s: (4) with localsearch=1 warns that it is ignored, should be one: %i \n',lab4,contains(warnmsg,'localsearch=1 is IGNORED'))
        fprintf('%s: (4) with localsearch=1 against without, V, this should be zero: %.3e \n',lab4,max(abs(VPath1(:)-VPath0(:))))
        fprintf('%s: (4) with localsearch=1 against without, Policy, this should be zero: %.3e \n',lab4,max(abs(PolicyPath1(:)-PolicyPath0(:))))
        vfoptionsF=vfoptionsL; vfoptionsF.localsearch=2;
        VPath2=ValueFnOnTransPath_Case1_FHorz(PricePathB, ParamPathGE, T, V_finalGE, Policy_finalGE, ParamsGE, n_d_cc, n_a_cc, n_z_cc, N_j, d_grid_cc, a_grid_cc,z_grid_cc, pi_z_cc, DiscountFactorParamNames, ReturnFn, transpathoptionsS, vfoptionsF);
        fprintf('%s: (4) with localsearch=2, nlocalsearch=1 (around staying put), V against without, should NOT be zero: %.3e \n',lab4,max(abs(VPath2(:)-VPath0(:))))
        fprintf('%s: (4) with localsearch=2, nlocalsearch=1, V is never above without, this should be zero: %.3e \n',lab4,max(max(VPath2(:)-VPath0(:)),0))
        vfoptionsF.nlocalsearch=(n_a_cc-1)/2;
        [VPath2c,PolicyPath2c]=ValueFnOnTransPath_Case1_FHorz(PricePathB, ParamPathGE, T, V_finalGE, Policy_finalGE, ParamsGE, n_d_cc, n_a_cc, n_z_cc, N_j, d_grid_cc, a_grid_cc,z_grid_cc, pi_z_cc, DiscountFactorParamNames, ReturnFn, transpathoptionsS, vfoptionsF);
        fprintf('%s: (4) with localsearch=2 at the cap, V against without, this should be zero: %.3e \n',lab4,max(abs(VPath2c(:)-VPath0(:))))
        fprintf('%s: (4) with localsearch=2 at the cap, Policy entries that differ (ties only): %i of %i \n',lab4,sum(PolicyPath2c(:)~=PolicyPath0(:)),numel(PolicyPath0))
    end
end

%% (5) Quasi-hyperbolic is still refused with local search (only Epstein-Zin was lifted)
vfoptionsE=struct(); vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=1;
vfoptionsE.exoticpreferences='QuasiHyperbolic'; vfoptionsE.QHadditionaldiscount='beta0';
ParamsE=ParamsGE; ParamsE.beta0=1; fired=0;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d_cc, n_a_cc, n_z_cc, N_j, d_grid_cc,a_grid_cc,z_grid_cc, pi_z_cc, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptionsE);
catch ME
    fired=contains(ME.message,'not implemented with exoticpreferences other than EpsteinZin');
end
fprintf('EZ localsearch FHorz: (5) local search with quasi-hyperbolic discounting is still refused by name, should be one: %i \n',fired)

output=1;

end
