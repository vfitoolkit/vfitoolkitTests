function output=EZFHorzTPath_d_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c)
% Epstein-Zin on the FHorz transition path, with d, with z, with e, without semiz.
% The EZ mirror of QHDFHorzTPath_d_z_e_nosemiz. Where QH runs the model twice (Naive, Sophisticated),
% this runs it three times, the three cases of CoreFHorzEZTests:
%   Case 1: consumption units (traditional Epstein-Zin, vfoptions.EZutils=0)
%   Case 2: utility units, positive-valued utility fn (EZutils=1, EZpositiveutility=1)
%   Case 3: utility units, negative-valued utility fn (EZutils=1, EZpositiveutility=0)
% For each case:
%   (A) the do-nothing path: prices and parameters held at their values in Params, so every period of the
%       path must equal the stationary Epstein-Zin solve, V and Policy, exactly. Plain, DC, GI and DC+GI,
%       both layouts (fastOLG=0 and 1), and every lowmemory the shocks allow. This pins each transition-
%       path raw to its stationary twin, which the stationary bank (CoreFHorzEZTests) has already pinned to
%       a brute-force oracle.
%   (B) on the moving path: fastOLG against slowOLG, DC against no DC, GI against DC+GI, lowmemory against not
%   (C) on the do-nothing path again, with the options that change the recursion: survival probabilities
%       with warm-glow of bequests; EZoneminusbeta=2 with EZmortalityriskaversion; EZoneminusbeta=1. Each
%       must be honoured in every period of the path, not dropped
%   (D) one iteration of TransitionPath_Case1_FHorz with a real general eqm condition: fastOLG against
%       slowOLG, and it must move the price (else the agreement is vacuous)
%   (F) the discriminator: on the moving path the Epstein-Zin V must differ from the vNM V of the same
%       return fn (should NOT be zero). Utility units only with shocks: without them the certainty
%       equivalent is the identity and the extra risk aversion has nothing to act on
% and once, after the three cases:
%   (E) the collapses on the moving path: consumption units with gamma=1/phi against vNM with the CRRA
%       return fn (Policy exact, V through the transform V_EZ=((1-gamma)*V_vNM).^(1/(1-gamma))); utility
%       units with EZriskaversion=0 against vNM with the same return fn (V and Policy exact); with and
%       without the grid interpolation layer
%   (G) case 1 on the big grid: moments with and without the grid interpolation layer (the figure)
% Until the transition path knew about Epstein-Zin, slowOLG solved an EZ model as if it were vNM, without
% a word, and fastOLG errored. (A) and (F) are the checks that fail on the first; slowOLG runs first in
% every loop so that they print before the second stops the run.

hasz=1; hase=1;
% with d: n_d and d_grid as passed in
% with z: n_z, z_grid and pi_z as passed in
vfoptionsS=struct();
simoptionsS=struct();
if hase==1
    vfoptionsS.n_e=vfoptionsbaseline.n_e;
    vfoptionsS.e_grid=vfoptionsbaseline.e_grid;
    vfoptionsS.pi_e=vfoptionsbaseline.pi_e;
    simoptionsS.n_e=simoptionsbaseline.n_e;
    simoptionsS.e_grid=simoptionsbaseline.e_grid;
    simoptionsS.pi_e=simoptionsbaseline.pi_e;
end
simoptionsSGI=simoptionsS; simoptionsSGI.gridinterplayer=1; simoptionsSGI.ngridinterp=5;

% zero assets, midpoints for any shocks
jsize=n_a; jsize_big=n_a_big; jmid={1};
if hasz==1
    jsize=[jsize,n_z]; jsize_big=[jsize_big,n_z]; jmid{end+1}=ceil(n_z/2);
end
if hase==1
    jsize=[jsize,vfoptionsS.n_e]; jsize_big=[jsize_big,vfoptionsS.n_e]; jmid{end+1}=ceil(vfoptionsS.n_e/2);
end
jequaloneDist=zeros([jsize,1],'gpuArray');
jequaloneDist(jmid{:})=1;
jequaloneDist_big=zeros([jsize_big,1],'gpuArray');
jequaloneDist_big(jmid{:})=1;

ReturnFn_cons=@(d,aprime,a,z,e,r,w,kappa_j,varphi,agej,Jr,pension) EZReturnFn_cons_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,varphi,agej,Jr,pension);
ReturnFn_posU=@(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension) EZReturnFn_positiveUtils_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
ReturnFn_negU=@(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension) EZReturnFn_negativeUtils_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
WarmGlowFn_cons=@(aprime,wg1,wg2) EZWarmGlowFn_cons(aprime,wg1,wg2);
WarmGlowFn_posU=@(aprime,wg1,wg2,wg3) EZWarmGlowFn_positiveUtils(aprime,wg1,wg2,wg3);
WarmGlowFn_negU=@(aprime,wg1,wg2,wg3) EZWarmGlowFn_negativeUtils(aprime,wg1,wg2,wg3);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(d,aprime,a,z,e) a;
FnsToEvaluate.earnings=@(d,aprime,a,z,e,w,kappa_j) w*kappa_j*z*e*d;

% A general eqm condition that depends on the agents' choices, for (D)
GeneralEqmEqns.CapitalMarket=@(r,assets) r-(0.3*((assets+0.1)^(-0.7))-0.08);

transpathoptionsL=cell(1,2); % slowOLG first, see the header
transpathoptionsL{1}=transpathoptionsbaseline; transpathoptionsL{1}.fastOLG=0;
transpathoptionsL{2}=transpathoptionsbaseline; transpathoptionsL{2}.fastOLG=1;
layoutname={'slowOLG','fastOLG'};

% Prices and parameters held at their values in Params: the do-nothing path
PricePathConstant.r=Params.r*ones(1,T);
ParamPathConstant.w=Params.w*ones(1,T);

% The lowmemory settings the shocks allow (lowmemory=2 loops over both z and e)
lowmemlist=[];
if hasz==1 || hase==1
    lowmemlist=1;
end
if hasz==1 && hase==1
    lowmemlist=[1,2];
end

tiername={'plain','DC','GI','DC+GI'};
casename={'cons units','positive utils','negative utils'};

for cc=1:3
    vfoptions=vfoptionsS;
    vfoptions.exoticpreferences='EpsteinZin';
    if cc==1
        vfoptions.EZutils=0; vfoptions.EZriskaversion='ezgamma'; vfoptions.EZeis='ezphi';
        ReturnFn=ReturnFn_cons; WarmGlowFn=WarmGlowFn_cons;
    elseif cc==2
        vfoptions.EZutils=1; vfoptions.EZpositiveutility=1; vfoptions.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU; WarmGlowFn=WarmGlowFn_posU;
    else
        vfoptions.EZutils=1; vfoptions.EZpositiveutility=0; vfoptions.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU; WarmGlowFn=WarmGlowFn_negU;
    end
    fprintf('\n--- EZ-TPath (d_z_e): %s ---\n',casename{cc})
    vfoptionsT=cell(1,4); % plain, DC, GI, DC+GI
    vfoptionsT{1}=vfoptions;
    vfoptionsT{2}=vfoptions; vfoptionsT{2}.divideandconquer=1;
    vfoptionsT{3}=vfoptions; vfoptionsT{3}.gridinterplayer=1; vfoptionsT{3}.ngridinterp=5;
    vfoptionsT{4}=vfoptions; vfoptionsT{4}.divideandconquer=1; vfoptionsT{4}.gridinterplayer=1; vfoptionsT{4}.ngridinterp=5;

    %% (A) The do-nothing path equals the stationary solve, period by period
    for tier=1:4
        [V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsT{tier});
        Vrep=reshape(repmat(V_s,[ones(1,ndims(V_s)),T]),[],1);
        Policyrep=reshape(repmat(Policy_s,[ones(1,ndims(Policy_s)),T]),[],1);
        for ll=1:2
            for lm=[0,lowmemlist]
                vfoptionsA=vfoptionsT{tier}; vfoptionsA.lowmemory=lm;
                [VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePathConstant, ParamPathConstant, T, V_s, Policy_s, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsL{ll}, vfoptionsA);
                fprintf('(A) %s: do-nothing path against the stationary solve (%s, %s, lowmemory=%i), this should be zero, V: %.3e \n',casename{cc},tiername{tier},layoutname{ll},lm,max(abs(VPath(:)-Vrep)))
                if numel(PolicyPath)==numel(Policyrep)
                    fprintf('(A) %s: do-nothing path against the stationary solve (%s, %s, lowmemory=%i), this should be zero, Policy: %.3e \n',casename{cc},tiername{tier},layoutname{ll},lm,max(abs(PolicyPath(:)-Policyrep)))
                else
                    fprintf('(A) %s: do-nothing path Policy is not the stationary Policy repeated over the path (%s, %s, lowmemory=%i), this should be zero, Policy: %.3e \n',casename{cc},tiername{tier},layoutname{ll},lm,1)
                end
            end
        end
    end
    clear VPath PolicyPath Vrep Policyrep

    %% (B) The moving path: layouts, divide-and-conquer, lowmemory
    [V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsT{1});
    [V_sGI,Policy_sGI]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsT{3});
    VfinalT={V_s,V_s,V_sGI,V_sGI};
    PolicyfinalT={Policy_s,Policy_s,Policy_sGI,Policy_sGI};
    VPathB=cell(2,4);
    PolicyPathB=cell(2,4);
    for ll=1:2
        for tier=1:4
            [VPathB{ll,tier},PolicyPathB{ll,tier}]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, VfinalT{tier}, PolicyfinalT{tier}, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsL{ll}, vfoptionsT{tier});
        end
    end
    for tier=1:4
        fprintf('(B) %s: fastOLG against slowOLG (%s), this should be zero, V: %.3e \n',casename{cc},tiername{tier},max(abs(VPathB{1,tier}(:)-VPathB{2,tier}(:))))
        fprintf('(B) %s: fastOLG against slowOLG (%s), this should be zero, Policy: %.3e \n',casename{cc},tiername{tier},max(abs(PolicyPathB{1,tier}(:)-PolicyPathB{2,tier}(:))))
    end
    for ll=1:2
        fprintf('(B) %s: divide-and-conquer against not (%s), this should be zero, V: %.3e \n',casename{cc},layoutname{ll},max(abs(VPathB{ll,1}(:)-VPathB{ll,2}(:))))
        fprintf('(B) %s: divide-and-conquer against not (%s), this should be zero, Policy: %.3e \n',casename{cc},layoutname{ll},max(abs(PolicyPathB{ll,1}(:)-PolicyPathB{ll,2}(:))))
        fprintf('(B) %s: divide-and-conquer against not (with GI, %s), this should be zero, V: %.3e \n',casename{cc},layoutname{ll},max(abs(VPathB{ll,3}(:)-VPathB{ll,4}(:))))
        fprintf('(B) %s: divide-and-conquer against not (with GI, %s), this should be zero, Policy: %.3e \n',casename{cc},layoutname{ll},max(abs(PolicyPathB{ll,3}(:)-PolicyPathB{ll,4}(:))))
    end
    % anti-vacuity: the moving path must move V off the stationary solve, or (B) compares two do-nothing paths
    Vrep=reshape(repmat(V_s,[ones(1,ndims(V_s)),T]),[],1);
    fprintf('(B) %s: the moving path moves V off the stationary solve, should NOT be zero: %.3e \n',casename{cc},max(abs(VPathB{1,1}(:)-Vrep)))
    clear Vrep
    for lm=lowmemlist
        for ll=1:2
            for tier=1:4
                vfoptionsB=vfoptionsT{tier}; vfoptionsB.lowmemory=lm;
                [VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, VfinalT{tier}, PolicyfinalT{tier}, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsL{ll}, vfoptionsB);
                fprintf('(B) %s: lowmemory=%i against lowmemory=0 (%s, %s), this should be zero, V: %.3e \n',casename{cc},lm,tiername{tier},layoutname{ll},max(abs(VPath(:)-VPathB{ll,tier}(:))))
                fprintf('(B) %s: lowmemory=%i against lowmemory=0 (%s, %s), this should be zero, Policy: %.3e \n',casename{cc},lm,tiername{tier},layoutname{ll},max(abs(PolicyPath(:)-PolicyPathB{ll,tier}(:))))
            end
        end
    end
    clear VPath PolicyPath

    %% (F) Discriminator: the Epstein-Zin path is not the vNM path of the same return fn
    if cc==1 || hasz==1 || hase==1
        vfoptionsN=vfoptionsS; % exoticpreferences not set: vNM
        [V_sN,Policy_sN]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsN);
        for ll=1:2
            VPathN=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_sN, Policy_sN, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsL{ll}, vfoptionsN);
            fprintf('(F) %s: Epstein-Zin path against the vNM path of the same return fn (%s), should NOT be zero: %.3e \n',casename{cc},layoutname{ll},max(abs(VPathB{ll,1}(:)-VPathN(:))))
        end
        clear VPathN V_sN Policy_sN
    end
    clear VPathB PolicyPathB

    %% (C) Options that change the recursion, on the do-nothing path
    % (C1) survival probabilities with warm-glow of bequests: plain and DC+GI, both layouts, and lowmemory
    vfoptionsC=vfoptions;
    vfoptionsC.survivalprobability='sj';
    vfoptionsC.WarmGlowBequestsFn=WarmGlowFn;
    vfoptionsCT={vfoptionsC,vfoptionsC};
    vfoptionsCT{2}.divideandconquer=1; vfoptionsCT{2}.gridinterplayer=1; vfoptionsCT{2}.ngridinterp=5;
    tiernameC={'plain','DC+GI'};
    for tier=1:2
        [V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsCT{tier});
        Vrep=reshape(repmat(V_s,[ones(1,ndims(V_s)),T]),[],1);
        Policyrep=reshape(repmat(Policy_s,[ones(1,ndims(Policy_s)),T]),[],1);
        for ll=1:2
            for lm=[0,lowmemlist]
                vfoptionsA=vfoptionsCT{tier}; vfoptionsA.lowmemory=lm;
                [VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePathConstant, ParamPathConstant, T, V_s, Policy_s, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsL{ll}, vfoptionsA);
                fprintf('(C1) %s: survival probabilities and warm-glow, do-nothing path against the stationary solve (%s, %s, lowmemory=%i), this should be zero, V: %.3e \n',casename{cc},tiernameC{tier},layoutname{ll},lm,max(abs(VPath(:)-Vrep)))
                fprintf('(C1) %s: survival probabilities and warm-glow, do-nothing path against the stationary solve (%s, %s, lowmemory=%i), this should be zero, Policy: %.3e \n',casename{cc},tiernameC{tier},layoutname{ll},lm,max(abs(PolicyPath(:)-Policyrep)))
            end
        end
        if tier==1
            % anti-vacuity: the options must change the stationary solve, or (C1) is (A) again
            [V_s0,~]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
            fprintf('(C1) %s: survival probabilities and warm-glow change the stationary V, should NOT be zero: %.3e \n',casename{cc},max(abs(V_s(:)-V_s0(:))))
            clear V_s0
        end
    end
    % (C2) EZoneminusbeta=2 (the (1-sj*beta) term) with a distinct mortality risk aversion; (C3) EZoneminusbeta=1
    vfoptionsC2=vfoptions;
    vfoptionsC2.survivalprobability='sj';
    vfoptionsC2.EZoneminusbeta=2;
    vfoptionsC2.EZmortalityriskaversion='ezmrisk';
    vfoptionsC3=vfoptions;
    vfoptionsC3.EZoneminusbeta=1;
    vfoptionsC3.gridinterplayer=1; vfoptionsC3.ngridinterp=5;
    vfoptionsCC={vfoptionsC2,vfoptionsC3};
    labelCC={'(C2) EZoneminusbeta=2 with EZmortalityriskaversion (plain)','(C3) EZoneminusbeta=1 (GI)'};
    for kk=1:2
        [V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsCC{kk});
        Vrep=reshape(repmat(V_s,[ones(1,ndims(V_s)),T]),[],1);
        Policyrep=reshape(repmat(Policy_s,[ones(1,ndims(Policy_s)),T]),[],1);
        for ll=1:2
            [VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePathConstant, ParamPathConstant, T, V_s, Policy_s, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsL{ll}, vfoptionsCC{kk});
            fprintf('%s %s: do-nothing path against the stationary solve (%s), this should be zero, V: %.3e \n',labelCC{kk},casename{cc},layoutname{ll},max(abs(VPath(:)-Vrep)))
            fprintf('%s %s: do-nothing path against the stationary solve (%s), this should be zero, Policy: %.3e \n',labelCC{kk},casename{cc},layoutname{ll},max(abs(PolicyPath(:)-Policyrep)))
        end
    end
    clear VPath PolicyPath Vrep Policyrep

    %% (D) One iteration of the general eqm path
    [V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_s,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsS);
    transpathoptionsD=transpathoptionsbaseline;
    transpathoptionsD.GEnewprice=3; % the shooting algorithm
    transpathoptionsD.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1};
    transpathoptionsD.maxiter=1;
    transpathoptionsD.verbose=0;
    transpathoptionsD.fastOLG=0;
    PricePathDslow=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_s, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsD, simoptionsS, vfoptions);
    transpathoptionsD.fastOLG=1;
    PricePathDfast=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_s, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsD, simoptionsS, vfoptions);
    fprintf('(D) %s: one iteration of the general eqm path, fastOLG against slowOLG, this should be zero: %.3e \n',casename{cc},max(abs(PricePathDslow.r-PricePathDfast.r)))
    fprintf('(D) %s: one iteration of the general eqm path moves r, should NOT be zero: %.3e \n',casename{cc},max(abs(PricePathDslow.r-PricePath.r)))
    clear AgentDist_initial

    %% (G) Case 1 only: moments on the big grid, with and without the grid interpolation layer
    if cc==1
        [V_sb,Policy_sb]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsT{2});
        [V_sbGI,Policy_sbGI]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsT{4});
        [~,PolicyPath2b]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_sb, Policy_sb, Params, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsT{2});
        AgentDist_initial_big=StationaryDist_FHorz_Case1(jequaloneDist_big,AgeWeightParamNames,Policy_sb,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptionsS);
        AgentDistPath2=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_big, jequaloneDist_big, PricePath, ParamPath, PolicyPath2b, AgeWeightParamNames,n_d,n_a_big,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptionsS);
        AggVarsPath2=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath2, PolicyPath2b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, transpathoptionsbaseline, simoptionsS);
        clear PolicyPath2b
        [~,PolicyPath4b]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_sbGI, Policy_sbGI, Params, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsT{4});
        AgentDistPath4=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_big, jequaloneDist_big, PricePath, ParamPath, PolicyPath4b, AgeWeightParamNames,n_d,n_a_big,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptionsSGI);
        AggVarsPath4=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath4, PolicyPath4b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, transpathoptionsbaseline, simoptionsSGI);
        clear PolicyPath4b
        fprintf('(G) %s: with/without grid interp, should get much the same moments (for big a_grid) \n',casename{cc})
        fprintf('(G) %s: AgentDist path with/without grid interp, this should be close to zero: %.3e \n',casename{cc},max(abs(AgentDistPath2(:)-AgentDistPath4(:))))
        fig=figure(figure_c);
        subplot(2,1,1); plot(1:1:T,AggVarsPath2.earnings.Mean, 1:1:T,AggVarsPath4.earnings.Mean)
        title('Earnings Mean (Epstein-Zin, consumption units)')
        legend('DC','DC+GI')
        subplot(2,1,2); plot(1:1:T,AggVarsPath2.assets.Mean, 1:1:T,AggVarsPath4.assets.Mean)
        title('Assets Mean (Epstein-Zin, consumption units)')
        legend('DC','DC+GI')
        clear V_sb V_sbGI Policy_sb Policy_sbGI AgentDist_initial_big AgentDistPath2 AgentDistPath4 AggVarsPath2 AggVarsPath4
    end
end

%% (E) The collapses to vNM, on the moving path
vfoptionsGI=vfoptionsS; vfoptionsGI.gridinterplayer=1; vfoptionsGI.ngridinterp=5;
vfoptionsE={vfoptionsS,vfoptionsGI};
tiernameE={'plain','GI'};
for gi=1:2
    % (E1) consumption units with gamma=1/phi is vNM with the CRRA return fn (x^(1-gamma))/(1-gamma), which
    % is the negative-utils return fn with ezsigma=gamma
    ParamsE=Params;
    ParamsE.ezphi=1/ParamsE.ezgamma;
    ParamsE.ezsigma=ParamsE.ezgamma;
    vfEZ=vfoptionsE{gi};
    vfEZ.exoticpreferences='EpsteinZin'; vfEZ.EZutils=0; vfEZ.EZriskaversion='ezgamma'; vfEZ.EZeis='ezphi';
    vfN=vfoptionsE{gi};
    [V_sEZ,Policy_sEZ]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,ParamsE,DiscountFactorParamNames,[],vfEZ);
    [V_sN,Policy_sN]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,ParamsE,DiscountFactorParamNames,[],vfN);
    for ll=1:2
        [VPathEZ,PolicyPathEZ]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_sEZ, Policy_sEZ, ParamsE, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn_cons, transpathoptionsL{ll}, vfEZ);
        [VPathN,PolicyPathN]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_sN, Policy_sN, ParamsE, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn_negU, transpathoptionsL{ll}, vfN);
        % The transform sends a dead state (-Inf, no feasible choice: zero assets at the ages with kappa_j=0)
        % to zero, so V is compared on the finite set and the dead sets are compared on their own
        VPathNtransformed=((1-ParamsE.ezgamma)*VPathN).^(1/(1-ParamsE.ezgamma));
        finEZ=isfinite(VPathEZ);
        fprintf('(E1) EZ gamma=1/phi collapse on the path (%s, %s), Policy, this should be zero: %.3e \n',tiernameE{gi},layoutname{ll},max(abs(PolicyPathEZ(:)-PolicyPathN(:))))
        fprintf('(E1) EZ gamma=1/phi collapse on the path (%s, %s), V after the transform on the finite set (relative), this should be zero: %.3e \n',tiernameE{gi},layoutname{ll},max(abs(VPathEZ(finEZ)-VPathNtransformed(finEZ)))/max(abs(VPathEZ(finEZ))))
        fprintf('(E1) EZ gamma=1/phi collapse on the path (%s, %s), the dead sets are the same, this should be one: %i \n',tiernameE{gi},layoutname{ll},isequal(gather(VPathEZ==-Inf),gather(VPathN==-Inf)))
    end
    clear VPathEZ PolicyPathEZ VPathN PolicyPathN VPathNtransformed finEZ
    % (E2) utility units with EZriskaversion=0 is vNM with the same return fn, both signs
    ParamsE=Params;
    ParamsE.ezrisk=0;
    ReturnFnE={ReturnFn_posU,ReturnFn_negU};
    for pp=1:2
        vfEZ=vfoptionsE{gi};
        vfEZ.exoticpreferences='EpsteinZin'; vfEZ.EZutils=1; vfEZ.EZpositiveutility=2-pp; vfEZ.EZriskaversion='ezrisk';
        [V_sEZ,Policy_sEZ]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFnE{pp},ParamsE,DiscountFactorParamNames,[],vfEZ);
        [V_sN,Policy_sN]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFnE{pp},ParamsE,DiscountFactorParamNames,[],vfN);
        for ll=1:2
            [VPathEZ,PolicyPathEZ]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_sEZ, Policy_sEZ, ParamsE, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFnE{pp}, transpathoptionsL{ll}, vfEZ);
            [VPathN,PolicyPathN]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_sN, Policy_sN, ParamsE, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFnE{pp}, transpathoptionsL{ll}, vfN);
            fprintf('(E2) EZ (%s) with EZriskaversion=0 on the path (%s, %s), should give zero, V: %.3e \n',casename{pp+1},tiernameE{gi},layoutname{ll},max(abs(VPathEZ(:)-VPathN(:))))
            fprintf('(E2) EZ (%s) with EZriskaversion=0 on the path (%s, %s), should give zero, Policy: %.3e \n',casename{pp+1},tiernameE{gi},layoutname{ll},max(abs(PolicyPathEZ(:)-PolicyPathN(:))))
        end
    end
    clear VPathEZ PolicyPathEZ VPathN PolicyPathN V_sEZ Policy_sEZ V_sN Policy_sN
end

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
