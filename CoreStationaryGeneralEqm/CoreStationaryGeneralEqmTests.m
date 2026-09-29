% Implement core Stationary General Eqm tests of the VFI Toolkit.
%
% Baseline model: Aiyagari-style incomplete markets with endogenous labor and
% a government. Three general eqm prices (r,Tr,tau_c) solve three general eqm
% eqns (CapitalMarket, GovBudget, ConsTax). The wage w is hardcoded from r via
% the firm FOC (not a price); the labor tax rate tau and spending G are fixed.
% See the setup file for details.
%
% For each horizon we:
%  (i)  solve with fminalgo=1, 5, 8, 4 and confirm all give the same answer
%        (fminalgo=4 is CMA-ES, stochastic/lower-accuracy, so a looser match)
%  (ii) re-solve using all three kinds of parameter constraint, each applied to
%        one of the prices: constrain0to1 on r, constrainpositive on Tr,
%        constrainAtoB on tau_c (and then all three at once); confirm each
%        reproduces the unconstrained answer. constrainpositive is done three
%        times: with constrainpositivemethod left unset, with it set explicitly
%        to 'softplus' (cparam=log(1+exp(uparam))), and with it set to 'log'
%        (uparam=log(cparam)). Since 2026-09-12 'softplus' IS the default, so the
%        first two are the same transform and must agree exactly; that pair is
%        what catches a GE solver whose own default has silently reverted. The
%        log run is the only one that puts a non-default constrainpositivemethod
%        through a GE solve, so the bank also checks that it DIFFERS from the
%        softplus run: an exact zero there means the option never reached the
%        solver, and every invariance check would stay green while covering
%        nothing. Both transforms send the real line to (0,infty), so both must
%        land on the same equilibrium; they differ in the step sizes the
%        optimizer sees, which matters when a constrained price starts at 1 (log
%        sends that to uparam=0, exactly where fminsearch switches to its tiny
%        absolute simplex perturbation). log is the worse of the two on this
%        model, so its invariance check carries a deliberately loose threshold
%        and the accuracy ratio is printed alongside it.
%
% This file runs: InfHorz, then FHorz, then the same again with permanent types
% (PType, N_i=2 types differing in sigma=2.2 and 1.8). Outputs: output1..output4
% for the non-PType tests, output1ptype..output4ptype for the PType ones.

% No figures are drawn anywhere in this bank, so only the diary is saved.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): constrainpositivemethod round trips
%            (pure transform checks - no model and no GE solve, so this part costs seconds;
%             the fminalgo and constraints parts are GE solves: 6 for each fminalgo part, 7 for
%             each constraints part, 52 in all. The extraoptions parts (10 to 13) are mostly NOT
%             solves: they use heteroagentoptions.maxiter=0, which evaluates the general eqm
%             conditions at the prices in Params and returns them, one model solve rather than a
%             hundred and fifty. Six real solves in the four of them, all told)
% doPart(2): InfHorz fminalgo agreement
% doPart(3): InfHorz parameter-constraint invariance
% doPart(4): FHorz fminalgo agreement
% doPart(5): FHorz parameter-constraint invariance
% doPart(6): InfHorz PType fminalgo agreement
% doPart(7): InfHorz PType parameter-constraint invariance
% doPart(8): FHorz PType fminalgo agreement
% doPart(9): FHorz PType parameter-constraint invariance
% doPart(10): InfHorz intermediateEqns and shock grids in GE
% doPart(11): FHorz intermediateEqns, jequaloneDist as a function, shock grids in GE
% doPart(12): InfHorz PType intermediateEqns and shock grids in GE
% doPart(13): FHorz PType intermediateEqns, the four jequaloneDist forms, shock grids in GE
% doPart(14): FHorz PType general eqm conditions by ptype (GEptype)
% doPart(15): InfHorz PType general eqm conditions by ptype (GEptype)
%             [ON again 2026-09-23: HeteroAgentStationaryEqm_InfHorz_PType_GEptype_subfn now exists,
%              so the 'not implemented for infinite horizon' error is gone. This part is its first run]
% doPart(16): FHorz PType GEptype parameter-constraint invariance (fminalgo 1, 5 and 9)
% doPart(17): InfHorz PType GEptype parameter-constraint invariance (fminalgo 1, 5 and 9)
%             [ON again 2026-09-23, same reason as doPart(15)]
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1];

%% Diary of the command window output
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreStationaryGeneralEqmTestsdiary.txt','file')
    delete('./TestOutput/CoreStationaryGeneralEqmTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreStationaryGeneralEqmTestsdiary.txt
fprintf('CoreStationaryGeneralEqmTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

%%
addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreStationaryGeneralEqm_subcodes/')
addpath('./CoreStationaryGeneralEqm_Setup/')
addpath('./CoreStationaryGeneralEqm_ReturnFns/')

% Setup so that we use the same model in both halves
CoreStationaryGE_setup


%% ===== doPart(1): constrainpositivemethod round trips =====
if doPart(1)==1
    fprintf('\n===== doPart(1): constrainpositivemethod round trips =====\n')
    %% ============ constrainpositivemethod round trips ============
    % Pure transform checks, no model and no GE solve, run first because if the
    % transforms are wrong there is no point solving anything with them.
    %
    % The GE-level checks further down put constrainpositive on Tr under both
    % methods, but Tr sits at a moderate value there, so those solves cannot see
    % what happens at the extremes. These three checks pin the parts that a GE solve
    % on this model would silently pass:
    %   (a) round trip over 16 orders of magnitude, both methods
    %   (b) cparam=100 under softplus. The +-50 cutoffs sit in log units, where +50
    %       means cparam of about 5e21; under softplus cparam is roughly uparam, so
    %       restoring an upper cutoff there would cap EVERY positive parameter at 50.
    %       This is the regression guard for that.
    %   (c) the overflow ends, where the naive log(exp(c)-1) and log(1+exp(u)) blow up.
    % Printed as %.3e because several of these differences are genuinely near a ULP.

    fprintf('\n=== constrainpositivemethod: transform round trips ===\n')
    cpm_names={'theta'};
    cpm_index=[0;1]; % one parameter, one element
    cpm_log.constrainpositive={'theta'}; cpm_log.constrain0to1={}; cpm_log.constrainAtoB={};
    cpm_log.constrainpositivemethod='log';
    cpm_sp=cpm_log; cpm_sp.constrainpositivemethod='softplus';

    % (a) round trip, and (b) is the 60/100/1e4 part of this grid
    cpm_grid=[1e-8,1e-4,0.01,0.5,1,2,10,50,60,100,1e4,1e8];
    cpm_worst_log=0; cpm_worst_sp=0;
    for cpm_c=cpm_grid
        [cpm_u1,cpm_o1]=ParameterConstraints_TransformParamsToUnconstrained(cpm_c,cpm_index,cpm_names,cpm_log,1);
        cpm_b1=ParameterConstraints_TransformParamsToOriginal(cpm_u1,cpm_index,cpm_names,cpm_o1);
        [cpm_u2,cpm_o2]=ParameterConstraints_TransformParamsToUnconstrained(cpm_c,cpm_index,cpm_names,cpm_sp,1);
        cpm_b2=ParameterConstraints_TransformParamsToOriginal(cpm_u2,cpm_index,cpm_names,cpm_o2);
        cpm_worst_log=max(cpm_worst_log,abs(cpm_b1-cpm_c)/cpm_c);
        cpm_worst_sp=max(cpm_worst_sp,abs(cpm_b2-cpm_c)/cpm_c);
        fprintf('cparam=%9.2e | log uparam=%10.4f back=%9.2e | softplus uparam=%10.4f back=%9.2e \n',cpm_c,cpm_u1,cpm_b1,cpm_u2,cpm_b2)
    end
    fprintf('worst round-trip relative error, log, this should be near zero: %.3e \n',cpm_worst_log)
    fprintf('worst round-trip relative error, softplus, this should be near zero: %.3e \n',cpm_worst_sp)
    fprintf('(the 60, 100 and 1e4 rows are the cutoff guard: if softplus ever regains \n')
    fprintf(' an upper cutoff they come back as about 50 instead of themselves) \n')

    % cparam=1 is why the option exists: log sends it to exactly uparam=0, which is
    % where fminsearch swaps its 5% relative simplex perturbation for an absolute
    % 0.00025, making the first step on that parameter about 200x too small.
    [cpm_u1,~]=ParameterConstraints_TransformParamsToUnconstrained(1,cpm_index,cpm_names,cpm_log,1);
    [cpm_u2,~]=ParameterConstraints_TransformParamsToUnconstrained(1,cpm_index,cpm_names,cpm_sp,1);
    fprintf('cparam=1: log uparam=%.6f (expect exactly 0), softplus uparam=%.6f (expect %.6f) \n',cpm_u1,cpm_u2,log(exp(1)-1))
    fprintf('log uparam(1) distance from 0, this should be zero: %.3e \n',abs(cpm_u1))
    fprintf('softplus uparam(1) distance from log(e-1), this should be near zero: %.3e \n',abs(cpm_u2-log(exp(1)-1)))

    % (c) overflow ends
    [cpm_ubig,cpm_o2]=ParameterConstraints_TransformParamsToUnconstrained(1e6,cpm_index,cpm_names,cpm_sp,1);
    cpm_bbig=ParameterConstraints_TransformParamsToOriginal(cpm_ubig,cpm_index,cpm_names,cpm_o2);
    cpm_chuge=ParameterConstraints_TransformParamsToOriginal(1e6,cpm_index,cpm_names,cpm_o2);
    fprintf('softplus forward cparam=1e6: uparam=%.6e back=%.6e relerr=%.3e \n',cpm_ubig,cpm_bbig,abs(cpm_bbig-1e6)/1e6)
    fprintf('softplus inverse uparam=1e6: cparam=%.6e \n',cpm_chuge)
    fprintf('count of Inf/NaN across the overflow checks, this should be zero: %d \n',sum(~isfinite([cpm_ubig,cpm_bbig,cpm_chuge])))

    % The default is 'softplus' (changed from 'log' on 2026-09-12, after the InfHorz
    % constraint-invariance check below came out at 5.3e-2 under log against 2.6e-5
    % under softplus; that comparison is no longer a one-off, the constraints
    % subcodes below now run both transforms through a GE solve and print the ratio
    % of their errors every time). So the option unset must equal 'softplus', and must NOT equal
    % 'log' -- the second half is what catches a default that has silently reverted.
    cpm_unset=rmfield(cpm_log,'constrainpositivemethod');
    cpm_maxdiff_sp=0; cpm_maxdiff_log=0;
    for cpm_c=cpm_grid
        [cpm_ua,~]=ParameterConstraints_TransformParamsToUnconstrained(cpm_c,cpm_index,cpm_names,cpm_unset,1);
        [cpm_ub,~]=ParameterConstraints_TransformParamsToUnconstrained(cpm_c,cpm_index,cpm_names,cpm_sp,1);
        [cpm_uc,~]=ParameterConstraints_TransformParamsToUnconstrained(cpm_c,cpm_index,cpm_names,cpm_log,1);
        cpm_maxdiff_sp=max(cpm_maxdiff_sp,abs(cpm_ua-cpm_ub));
        cpm_maxdiff_log=max(cpm_maxdiff_log,abs(cpm_ua-cpm_uc));
    end
    fprintf('option unset vs method=''softplus'', this should be zero: %.3e \n',cpm_maxdiff_sp)
    fprintf('option unset vs method=''log'', this should NOT be zero: %.3e \n',cpm_maxdiff_log)
end % doPart(1): constrainpositivemethod round trips


%% ========================= InfHorz ==========================

%% ===== doPart(2): InfHorz fminalgo agreement =====
if doPart(2)==1
    fprintf('\n===== doPart(2): InfHorz fminalgo agreement =====\n')
    %% (i)  fminalgo agreement
    output1=CoreStationaryGE_InfHorz_fminalgo(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);

    % fprintf('\n=== InfHorz: fminalgo agreement (r,Tr,tau_c) ===\n')
    % fprintf('fminalgo=1: r=%.6f Tr=%.6f tau_c=%.6f \n',output1.p_eqm1.r,output1.p_eqm1.Tr,output1.p_eqm1.tau_c)
    % fprintf('fminalgo=5: r=%.6f Tr=%.6f tau_c=%.6f \n',output1.p_eqm5.r,output1.p_eqm5.Tr,output1.p_eqm5.tau_c)
    % fprintf('fminalgo=8: r=%.6f Tr=%.6f tau_c=%.6f \n',output1.p_eqm8.r,output1.p_eqm8.Tr,output1.p_eqm8.tau_c)
    % fprintf('fminalgo=4: r=%.6f Tr=%.6f tau_c=%.6f \n',output1.p_eqm4.r,output1.p_eqm4.Tr,output1.p_eqm4.tau_c)
    % d15=max(abs([output1.p_eqm1.r-output1.p_eqm5.r,output1.p_eqm1.Tr-output1.p_eqm5.Tr,output1.p_eqm1.tau_c-output1.p_eqm5.tau_c]));
    % d18=max(abs([output1.p_eqm1.r-output1.p_eqm8.r,output1.p_eqm1.Tr-output1.p_eqm8.Tr,output1.p_eqm1.tau_c-output1.p_eqm8.tau_c]));
    % d14=max(abs([output1.p_eqm1.r-output1.p_eqm4.r,output1.p_eqm1.Tr-output1.p_eqm4.Tr,output1.p_eqm1.tau_c-output1.p_eqm4.tau_c]));
    % fprintf('fminalgo 1 vs 5, this should be near zero: %.8f \n',d15)
    % fprintf('fminalgo 1 vs 8, this should be near zero: %.8f \n',d18)
    % fprintf('fminalgo 1 vs 4, this should be near zero: %.8f \n',d14)
end % doPart(2): InfHorz fminalgo agreement


%% ===== doPart(3): InfHorz parameter-constraint invariance =====
if doPart(3)==1
    fprintf('\n===== doPart(3): InfHorz parameter-constraint invariance =====\n')
    %% (ii) parameter-constraint invariance
    output2=CoreStationaryGE_InfHorz_constraints(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);

    % fprintf('\n=== InfHorz: parameter-constraint invariance ===\n')
    % fprintf('unconstrained:        r=%.6f Tr=%.6f tau_c=%.6f \n',output2.p_eqm0.r,output2.p_eqm0.Tr,output2.p_eqm0.tau_c)
    % fprintf('constrain0to1 on r:   r=%.6f Tr=%.6f tau_c=%.6f \n',output2.p_eqmA.r,output2.p_eqmA.Tr,output2.p_eqmA.tau_c)
    % fprintf('constrainpositive Tr: r=%.6f Tr=%.6f tau_c=%.6f \n',output2.p_eqmB.r,output2.p_eqmB.Tr,output2.p_eqmB.tau_c)
    % fprintf('constrainAtoB tau_c:  r=%.6f Tr=%.6f tau_c=%.6f \n',output2.p_eqmC.r,output2.p_eqmC.Tr,output2.p_eqmC.tau_c)
    % fprintf('constrain all three:  r=%.6f Tr=%.6f tau_c=%.6f \n',output2.p_eqmD.r,output2.p_eqmD.Tr,output2.p_eqmD.tau_c)
    % dA=max(abs([output2.p_eqm0.r-output2.p_eqmA.r,output2.p_eqm0.Tr-output2.p_eqmA.Tr,output2.p_eqm0.tau_c-output2.p_eqmA.tau_c]));
    % dB=max(abs([output2.p_eqm0.r-output2.p_eqmB.r,output2.p_eqm0.Tr-output2.p_eqmB.Tr,output2.p_eqm0.tau_c-output2.p_eqmB.tau_c]));
    % dC=max(abs([output2.p_eqm0.r-output2.p_eqmC.r,output2.p_eqm0.Tr-output2.p_eqmC.Tr,output2.p_eqm0.tau_c-output2.p_eqmC.tau_c]));
    % dD=max(abs([output2.p_eqm0.r-output2.p_eqmD.r,output2.p_eqm0.Tr-output2.p_eqmD.Tr,output2.p_eqm0.tau_c-output2.p_eqmD.tau_c]));
    % fprintf('constrain0to1 on r, this should be near zero: %.8f \n',dA)
    % fprintf('constrainpositive on Tr, this should be near zero: %.8f \n',dB)
    % fprintf('constrainAtoB on tau_c, this should be near zero: %.8f \n',dC)
    % fprintf('constrain all three, this should be near zero: %.8f \n',dD)
end % doPart(3): InfHorz parameter-constraint invariance





%% ========================== FHorz ===========================

%% ===== doPart(4): FHorz fminalgo agreement =====
if doPart(4)==1
    fprintf('\n===== doPart(4): FHorz fminalgo agreement =====\n')
    %% (i)  fminalgo agreement
    output3=CoreStationaryGE_FHorz_fminalgo(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);

    % fprintf('\n=== InfHorz: fminalgo agreement (r,Tr,tau_c) ===\n')
    % fprintf('fminalgo=1: r=%.6f Tr=%.6f tau_c=%.6f \n',output3.p_eqm1.r,output3.p_eqm1.Tr,output3.p_eqm1.tau_c)
    % fprintf('fminalgo=5: r=%.6f Tr=%.6f tau_c=%.6f \n',output3.p_eqm5.r,output3.p_eqm5.Tr,output3.p_eqm5.tau_c)
    % fprintf('fminalgo=8: r=%.6f Tr=%.6f tau_c=%.6f \n',output3.p_eqm8.r,output3.p_eqm8.Tr,output3.p_eqm8.tau_c)
    % fprintf('fminalgo=4: r=%.6f Tr=%.6f tau_c=%.6f \n',output3.p_eqm4.r,output3.p_eqm4.Tr,output3.p_eqm4.tau_c)
    % d15=max(abs([output3.p_eqm1.r-output3.p_eqm5.r,output3.p_eqm1.Tr-output3.p_eqm5.Tr,output3.p_eqm1.tau_c-output3.p_eqm5.tau_c]));
    % d18=max(abs([output3.p_eqm1.r-output3.p_eqm8.r,output3.p_eqm1.Tr-output3.p_eqm8.Tr,output3.p_eqm1.tau_c-output3.p_eqm8.tau_c]));
    % d14=max(abs([output3.p_eqm1.r-output3.p_eqm4.r,output3.p_eqm1.Tr-output3.p_eqm4.Tr,output3.p_eqm1.tau_c-output3.p_eqm4.tau_c]));
    % fprintf('fminalgo 1 vs 5, this should be near zero: %.8f \n',d15)
    % fprintf('fminalgo 1 vs 8, this should be near zero: %.8f \n',d18)
    % fprintf('fminalgo 1 vs 4, this should be near zero: %.8f \n',d14)
end % doPart(4): FHorz fminalgo agreement

%% ===== doPart(5): FHorz parameter-constraint invariance =====
if doPart(5)==1
    fprintf('\n===== doPart(5): FHorz parameter-constraint invariance =====\n')
    %% (ii) parameter-constraint invariance
    output4=CoreStationaryGE_FHorz_constraints(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);

    % fprintf('\n=== FHorz: parameter-constraint invariance ===\n')
    % fprintf('unconstrained:        r=%.6f Tr=%.6f tau_c=%.6f \n',output4.p_eqm0.r,output4.p_eqm0.Tr,output4.p_eqm0.tau_c)
    % fprintf('constrain0to1 on r:   r=%.6f Tr=%.6f tau_c=%.6f \n',output4.p_eqmA.r,output4.p_eqmA.Tr,output4.p_eqmA.tau_c)
    % fprintf('constrainpositive Tr: r=%.6f Tr=%.6f tau_c=%.6f \n',output4.p_eqmB.r,output4.p_eqmB.Tr,output4.p_eqmB.tau_c)
    % fprintf('constrainAtoB tau_c:  r=%.6f Tr=%.6f tau_c=%.6f \n',output4.p_eqmC.r,output4.p_eqmC.Tr,output4.p_eqmC.tau_c)
    % fprintf('constrain all three:  r=%.6f Tr=%.6f tau_c=%.6f \n',output4.p_eqmD.r,output4.p_eqmD.Tr,output4.p_eqmD.tau_c)
    % dA=max(abs([output4.p_eqm0.r-output4.p_eqmA.r,output4.p_eqm0.Tr-output4.p_eqmA.Tr,output4.p_eqm0.tau_c-output4.p_eqmA.tau_c]));
    % dB=max(abs([output4.p_eqm0.r-output4.p_eqmB.r,output4.p_eqm0.Tr-output4.p_eqmB.Tr,output4.p_eqm0.tau_c-output4.p_eqmB.tau_c]));
    % dC=max(abs([output4.p_eqm0.r-output4.p_eqmC.r,output4.p_eqm0.Tr-output4.p_eqmC.Tr,output4.p_eqm0.tau_c-output4.p_eqmC.tau_c]));
    % dD=max(abs([output4.p_eqm0.r-output4.p_eqmD.r,output4.p_eqm0.Tr-output4.p_eqmD.Tr,output4.p_eqm0.tau_c-output4.p_eqmD.tau_c]));
    % fprintf('constrain0to1 on r, this should be near zero: %.8f \n',dA)
    % fprintf('constrainpositive on Tr, this should be near zero: %.8f \n',dB)
    % fprintf('constrainAtoB on tau_c, this should be near zero: %.8f \n',dC)
    % fprintf('constrain all three, this should be near zero: %.8f \n',dD)
end % doPart(5): FHorz parameter-constraint invariance




%% ===================== InfHorz with PType ===================
% Same tests again, but with N_i=2 permanent types differing in sigma (2.2 and 1.8).

%% ===== doPart(6): InfHorz PType fminalgo agreement =====
if doPart(6)==1
    fprintf('\n===== doPart(6): InfHorz PType fminalgo agreement =====\n')
    %% (i)  fminalgo agreement
    output1ptype=CoreStationaryGE_InfHorz_PType_fminalgo(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);
end % doPart(6): InfHorz PType fminalgo agreement

%% ===== doPart(7): InfHorz PType parameter-constraint invariance =====
if doPart(7)==1
    fprintf('\n===== doPart(7): InfHorz PType parameter-constraint invariance =====\n')
    %% (ii) parameter-constraint invariance
    output2ptype=CoreStationaryGE_InfHorz_PType_constraints(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);
end % doPart(7): InfHorz PType parameter-constraint invariance





%% ====================== FHorz with PType ====================
% Same tests again, but with N_i=2 permanent types differing in sigma (2.2 and 1.8).

%% ===== doPart(8): FHorz PType fminalgo agreement =====
if doPart(8)==1
    fprintf('\n===== doPart(8): FHorz PType fminalgo agreement =====\n')
    %% (i)  fminalgo agreement
    output3ptype=CoreStationaryGE_FHorz_PType_fminalgo(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);
end % doPart(8): FHorz PType fminalgo agreement

%% ===== doPart(9): FHorz PType parameter-constraint invariance =====
if doPart(9)==1
    fprintf('\n===== doPart(9): FHorz PType parameter-constraint invariance =====\n')
    %% (ii) parameter-constraint invariance
    output4ptype=CoreStationaryGE_FHorz_PType_constraints(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);
end % doPart(9): FHorz PType parameter-constraint invariance



% The 'redo all these tests with GEbyptype' note that used to sit here is done: doPart(14) and (15)
% cover general eqm conditions by ptype, and doPart(16) and (17) sweep the parameter constraints
% under them, under the three solvers that are worth distinguishing there (see those subcodes).
% The deliberate remaining gap is a full fminalgo sweep under GEptype: fminalgo 4 and 8 reach the
% equilibrium through the same objective function as 1, so they would only re-test that.

%% ===== doPart(10): InfHorz intermediateEqns and shock grids in GE =====
if doPart(10)==1
    fprintf('\n===== doPart(10): InfHorz intermediateEqns and shock grids in GE =====\n')
    output5=CoreStationaryGE_InfHorz_extraoptions(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);
end % doPart(10): InfHorz intermediateEqns and shock grids in GE

%% ===== doPart(11): FHorz intermediateEqns, jequaloneDist as a function, shock grids in GE =====
if doPart(11)==1
    fprintf('\n===== doPart(11): FHorz intermediateEqns, jequaloneDist as a function, shock grids in GE =====\n')
    output6=CoreStationaryGE_FHorz_extraoptions(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);
end % doPart(11): FHorz intermediateEqns, jequaloneDist as a function, shock grids in GE

%% ===== doPart(12): InfHorz PType intermediateEqns and shock grids in GE =====
if doPart(12)==1
    fprintf('\n===== doPart(12): InfHorz PType intermediateEqns and shock grids in GE =====\n')
    output5ptype=CoreStationaryGE_InfHorz_PType_extraoptions(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);
end % doPart(12): InfHorz PType intermediateEqns and shock grids in GE

%% ===== doPart(13): FHorz PType intermediateEqns, the four jequaloneDist forms, shock grids in GE =====
if doPart(13)==1
    fprintf('\n===== doPart(13): FHorz PType intermediateEqns, the four jequaloneDist forms, shock grids in GE =====\n')
    output6ptype=CoreStationaryGE_FHorz_PType_extraoptions(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);
end % doPart(13): FHorz PType intermediateEqns, the four jequaloneDist forms, shock grids in GE


%% ===== doPart(14): FHorz PType general eqm conditions by ptype (GEptype) =====
if doPart(14)==1
    fprintf('\n===== doPart(14): FHorz PType general eqm conditions by ptype (GEptype) =====\n')
    output7ptype=CoreStationaryGE_FHorz_PType_GEptype(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);
end % doPart(14): FHorz PType general eqm conditions by ptype (GEptype)

%% ===== doPart(15): InfHorz PType general eqm conditions by ptype (GEptype) =====
if doPart(15)==1
    fprintf('\n===== doPart(15): InfHorz PType general eqm conditions by ptype (GEptype) =====\n')
    output8ptype=CoreStationaryGE_InfHorz_PType_GEptype(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);
end % doPart(15): InfHorz PType general eqm conditions by ptype (GEptype)


%% ===== doPart(16): FHorz PType GEptype parameter-constraint invariance =====
if doPart(16)==1
    fprintf('\n===== doPart(16): FHorz PType GEptype parameter-constraint invariance =====\n')
    output9ptype=CoreStationaryGE_FHorz_PType_GEptype_constraints(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);
end % doPart(16): FHorz PType GEptype parameter-constraint invariance

%% ===== doPart(17): InfHorz PType GEptype parameter-constraint invariance =====
if doPart(17)==1
    fprintf('\n===== doPart(17): InfHorz PType GEptype parameter-constraint invariance =====\n')
    output10ptype=CoreStationaryGE_InfHorz_PType_GEptype_constraints(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);
end % doPart(17): InfHorz PType GEptype parameter-constraint invariance


%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreStationaryGeneralEqmTestsdiary.txt')

diary off
