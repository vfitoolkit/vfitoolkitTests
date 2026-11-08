function output=EZFHorzTPath_reform(T,n_d,n_a,n_z,n_e,N_j,d_grid,a_grid,z_grid,pi_z,e_grid,pi_e,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline)
% An actual reform with Epstein-Zin preferences: a permanent rise in risk aversion, announced at t=1. No figures.
% Every other general eqm path in this bank is a null reform (prices and parameters at their stationary values,
% solved from a bumped guess), so none of them shows that an Epstein-Zin parameter moving on the ParamPath gives
% the right general eqm. Here the economy starts in the stationary general eqm at the old risk aversion and must
% end in the stationary general eqm at the new one. The model is the one in EZFHorzTPath_localsearch_shooting
% (doPart(7)): a Cobb-Douglas firm supplies both r and w; survival probabilities and the warm-glow are on.
%
% Risk aversion is vfoptions.EZriskaversion: ezgamma for consumption units, ezrisk for the utility-units cases.
% It goes from its value in Params to that plus one, on the ParamPath in every period (and in Params, which is
% the final stationary eqm's).
%
% Per row of the configuration table below:
%   (0) the old and the new stationary general eqm are each within the path's toleranceGEcondns, and the
%       reform moves them (or the path would be a null reform again)
%   (1) the plain shooting path reaches toleranceGEcondns
%   (2) the path ends at the new stationary general eqm: its last periods' prices against the new
%       stationary prices; and it starts away from them (the old distribution is not the new one)
%   (3) the same path with local search: within tolerance, the same price path, and an INDEPENDENT plain sweep
%       at its prices is within tolerance too

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

% columns: EZ case (1=cons units, 2=positive utils, 3=negative utils), fastOLG, shocks (1=z,2=noz,3=e,4=noz_e), hasd, GI
config=[1,0,1,0,0;   % the three EZ cases, slowOLG, z
        2,1,1,0,1;
        3,0,1,0,0;
        3,1,4,0,0];  % e without z
casename={'cons units','positive utils','negative utils'};
casefn={'cons','positiveUtils','negativeUtils'};
layoutname={'slowOLG','fastOLG'};
shockname={'z','noz','e','noz_e'};

for cc=1:size(config,1)
    ezcase=config(cc,1); fastOLG=config(cc,2); ss=config(cc,3); hasd=config(cc,4); useGI=config(cc,5);
    hasz=(ss==1 || ss==3);
    hase=(ss==3 || ss==4);
    lab=sprintf('EZ reform FHorz (%s, %s, %s, d%i, GI%i)',casename{ezcase},layoutname{fastOLG+1},shockname{ss},hasd,useGI);

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
    n_a_cc=n_a; a_grid_cc=a_grid;
    if cc==1 || cc==3
        % slowOLG, z, without GI: at n_a=101 the cons units path stalled at 1.29e-4 after 1000 iterations (the
        % reform takes r from 0.287 to 0.250, prices the null reform never visits), so a finer grid of the same shape
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

    %% The old and the new stationary general eqm
    % As in doPart(7): solved to 1e-8, so that each sits inside the path's toleranceGEcondns (1e-4)
    heteroagentoptionsGE=struct();
    heteroagentoptionsGE.verbose=0;
    heteroagentoptionsGE.toleranceGEprices=1e-8;
    heteroagentoptionsGE.toleranceGEcondns=1e-8;
    heteroagentoptionsGE.maxiter=2000;
    heteroagentoptionsGE.outputGEstruct=2;
    % old: risk aversion as in Params
    ParamsOld=Params;
    [p_old,GEcondnsOld]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d_cc,n_a_cc,n_z_cc,N_j,0,pi_z_cc,d_grid_cc,a_grid_cc,z_grid_cc,ReturnFn,FnsToEvaluateGE,GeneralEqmEqnsGE,ParamsOld,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptionsGE,simoptions,vfoptions);
    ParamsOld.r=p_old.r;
    ParamsOld.w=p_old.w;
    [~,Policy_old]=ValueFnIter_Case1_FHorz(n_d_cc,n_a_cc,n_z_cc,N_j,d_grid_cc,a_grid_cc,z_grid_cc,pi_z_cc,ReturnFn,ParamsOld,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_old,n_d_cc,n_a_cc,n_z_cc,N_j,pi_z_cc,ParamsOld,simoptions);
    % new: risk aversion one higher, from the old prices as the initial guess
    ParamsNew=ParamsOld;
    ParamsNew.(vfoptions.EZriskaversion)=Params.(vfoptions.EZriskaversion)+1;
    [p_new,GEcondnsNew]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d_cc,n_a_cc,n_z_cc,N_j,0,pi_z_cc,d_grid_cc,a_grid_cc,z_grid_cc,ReturnFn,FnsToEvaluateGE,GeneralEqmEqnsGE,ParamsNew,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptionsGE,simoptions,vfoptions);
    ParamsNew.r=p_new.r;
    ParamsNew.w=p_new.w;
    [V_final,~]=ValueFnIter_Case1_FHorz(n_d_cc,n_a_cc,n_z_cc,N_j,d_grid_cc,a_grid_cc,z_grid_cc,pi_z_cc,ReturnFn,ParamsNew,DiscountFactorParamNames,[],vfoptions);
    % (0) root sum of squares, the same distance the path uses per period under multiGEcriterion=1
    distOld=sqrt(sum(gather(GEcondnsOld(:)).^2));
    distNew=sqrt(sum(gather(GEcondnsNew(:)).^2));
    fprintf('%s: (0) old stationary general eqm conditions, root sum of squares, within the path tolerance, should be one: %i (it is %8.3e) \n',lab,distOld<=1e-4,distOld) % 1e-4 is the path's toleranceGEcondns, set below
    fprintf('%s: (0) new stationary general eqm conditions, root sum of squares, within the path tolerance, should be one: %i (it is %8.3e) \n',lab,distNew<=1e-4,distNew)
    fprintf('%s: (0) the reform moves the stationary eqm, r new against old, should NOT be zero: %.3e (r old %.6f, new %.6f) \n',lab,abs(p_new.r-p_old.r),p_old.r,p_new.r)

    %% The path: risk aversion at its new value in every period, from the old distribution
    ParamPath=struct();
    ParamPath.(vfoptions.EZriskaversion)=ParamsNew.(vfoptions.EZriskaversion)*ones(1,T);
    % initial guess: from the old to the new prices over the first half, then the new
    Thalf=floor(T/2);
    PricePath0.r=[linspace(p_old.r,p_new.r,Thalf),p_new.r*ones(1,T-Thalf)];
    PricePath0.w=[linspace(p_old.w,p_new.w,Thalf),p_new.w*ones(1,T-Thalf)];

    transpathoptionsS=transpathoptionsbaseline;
    transpathoptionsS.fastOLG=fastOLG;
    transpathoptionsS.GEnewprice=3;
    transpathoptionsS.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
    transpathoptionsS.maxiter=1000;
    transpathoptionsS.verbose=0;
    % Set here, not left to the default, because the checks read them: the distance below is the one the
    % solver uses under multiGEcriterion=1 (root sum of squares per period, then max over periods)
    transpathoptionsS.toleranceGEcondns=1e-4;
    transpathoptionsS.multiGEcriterion=1;
    if fastOLG==1
        simoptions.fastOLG=1;
    end

    %% (1) the plain solve
    tic;
    [PricePathB,GEcondnPathB]=TransitionPath_Case1_FHorz(PricePath0, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d_cc, n_a_cc, n_z_cc, N_j, d_grid_cc,a_grid_cc,z_grid_cc, pi_z_cc, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsNew, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptions);
    timeB=toc;
    distB=max(sqrt(GEcondnPathB.CapitalMarket.^2+GEcondnPathB.LabourMarket.^2));
    fprintf('%s: (1) plain general eqm conditions within tolerance, should be one: %i \n',lab,distB<=transpathoptionsS.toleranceGEcondns)

    %% (2) where the path ends, and where it starts
    % The last periods are back at the new stationary general eqm, to about the general eqm tolerance
    Ttail=T-9:T;
    fprintf('%s: (2) the last ten periods against the new stationary eqm, r, should be close to zero: %.3e \n',lab,max(abs(PricePathB.r(Ttail)-p_new.r)))
    fprintf('%s: (2) the last ten periods against the new stationary eqm, w, should be close to zero: %.3e \n',lab,max(abs(PricePathB.w(Ttail)-p_new.w)))
    fprintf('%s: (2) the first period is away from the new stationary eqm, r, should NOT be zero: %.3e \n',lab,abs(PricePathB.r(1)-p_new.r))

    %% (3) the same solve with local search, certified
    vfoptionsL=vfoptions;
    vfoptionsL.localsearch=1;
    vfoptionsL.nlocalsearch=1;
    tic;
    [PricePathL,GEcondnPathL]=TransitionPath_Case1_FHorz(PricePath0, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d_cc, n_a_cc, n_z_cc, N_j, d_grid_cc,a_grid_cc,z_grid_cc, pi_z_cc, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsNew, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptionsL);
    timeL=toc;
    fprintf('%s: plain took %.2f seconds, with local search %.2f seconds \n',lab,timeB,timeL)
    distL=max(sqrt(GEcondnPathL.CapitalMarket.^2+GEcondnPathL.LabourMarket.^2));
    fprintf('%s: (3) local search general eqm conditions within tolerance, should be one: %i \n',lab,distL<=transpathoptionsS.toleranceGEcondns)
    % Both are solved to toleranceGEcondns from the same start, so they agree to about that tolerance
    fprintf('%s: (3) price path with local search against without, r, should be close to zero: %.3e \n',lab,max(abs(PricePathL.r-PricePathB.r)))
    fprintf('%s: (3) price path with local search against without, w, should be close to zero: %.3e \n',lab,max(abs(PricePathL.w-PricePathB.w)))
    % Independently of the solver's bookkeeping: one plain shooting iteration AT the local search answer
    transpathoptionsX=transpathoptionsS;
    transpathoptionsX.maxiter=1;
    [~,GEcondnPathX]=TransitionPath_Case1_FHorz(PricePathL, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d_cc, n_a_cc, n_z_cc, N_j, d_grid_cc,a_grid_cc,z_grid_cc, pi_z_cc, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsNew, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsX, simoptions, vfoptions);
    distX=max(sqrt(GEcondnPathX.CapitalMarket.^2+GEcondnPathX.LabourMarket.^2));
    fprintf('%s: (3) the local search answer re-evaluated by an independent plain sweep, general eqm conditions within tolerance, should be one: %i \n',lab,distX<=transpathoptionsS.toleranceGEcondns)
    % Shooting returns the conditions of its final (standard) sweep at the returned prices, so these agree exactly
    fprintf('%s: (3) local search conditions against the independent plain sweep at the same prices, this should be zero: %.3e \n',lab,max(abs(GEcondnPathL.CapitalMarket-GEcondnPathX.CapitalMarket))+max(abs(GEcondnPathL.LabourMarket-GEcondnPathX.LabourMarket)))
end

output=1;

end
