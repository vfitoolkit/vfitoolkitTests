% Implement lots of tests of the core VFI Toolkit InfHorz commands
% with/without d
% with/without z
% with/without grid interpolation
% with/without low memory (where appropriate; z cases only)
%
% NOT tested here:
%   - divide-and-conquer: not usable for InfHorz (too slow); it is tested in
%     the InfHorz-TPath test bank instead
%   - e/semiz: InfHorz core is just with/without d and with/without z here


%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): one endogenous state (figs 1-6)
% doPart(2): cross-tests
% doPart(3): two endogenous states (figs 17-18)
% doPart(4): two-endo cross-tests: second state ignored
% doPart(5): two-endo cross-tests: first state ignored
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreInfHorzTestsdiary.txt','file')
    delete('./TestOutput/CoreInfHorzTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreInfHorzTestsdiary.txt
fprintf('CoreInfHorzTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))


addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreInfHorzTests_subcodes/')
addpath('./CoreInfHorzTests_Setup/')
addpath('./CoreInfHorz_ReturnFns/')
addpath('./CoreInfHorzTests_subcodes/CrossTests/')
% Setup so that use the same d,a,z in all the models that use them
CoreInfHorz_setup

%% ===== doPart(1): one endogenous state (figs 1-6) =====
if doPart(1)==1
    fprintf('\n===== doPart(1): one endogenous state (figs 1-6) =====\n')
    %% without d, without z
    figure_c=1;
    output=CoreInfHorz_nod_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % Note: no exogenous shock, so the stationary dist is a single mass point and the distribution statistics are degenerate. 
    % Solver/GI/ValueFnFromPolicy tests are the meaningful ones here.

    %% with d, without z
    figure_c=2;
    output=CoreInfHorz_d_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, with z
    figure_c=3;
    output=CoreInfHorz_nod_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z
    figure_c=4;
    output=CoreInfHorz_d_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, with e (iid), without z
    figure_c=5;
    % output=CoreInfHorz_nod_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with e (iid), without z
    figure_c=6;
    % output=CoreInfHorz_d_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(1): one endogenous state (figs 1-6)

%% ===== doPart(2): cross-tests =====
if doPart(2)==1
    fprintf('\n===== doPart(2): cross-tests =====\n')
    %% Cross-tests: run things that should give the same answer, disguised different ways, and check they do
    % A: a single trivial z point (value 1, prob 1) reproduces the noz code path
    % B: an iid disguised as a markov (identical rows) => stationary marginal over z equals that row
    % C: relabelling (permuting) the z states leaves all aggregate statistics unchanged
    output=CoreInfHorz_CrossTests_nod_nosemiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreInfHorz_CrossTests_d_nosemiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(2): cross-tests

% % % %% Cross-tests for the iid e (no z) models
% % % % D: an iid e, disguised as a markov z with identical rows, gives the same answer
% % % % E: a single trivial e point (value 1, prob 1) reproduces the noe code path
% % % % F: relabelling (permuting) the e states leaves all aggregate statistics unchanged
% % % output=CoreInfHorz_CrossTests_nod_nosemiz_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
% % % 
% % % output=CoreInfHorz_CrossTests_d_nosemiz_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);

%% Two endogenous states (mirror figs 3 & 4, but with a second endogenous state)
% Model: fig 17 simplified Kitao (2008), no taxes -- assets + occupation (worker/entrepreneur).
% fig 18 adds Bruggemann (2021) endogenous labor supply as the decision variable d.

%% ===== doPart(3): two endogenous states (figs 17-18) =====
if doPart(3)==1
    fprintf('\n===== doPart(3): two endogenous states (figs 17-18) =====\n')
    %% without d, with z, two endogenous states
    figure_c=17;
    output=CoreInfHorz_nod_z_noe_nosemiz_with2A(Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d (endogenous labor), with z, two endogenous states
    figure_c=18;
    output=CoreInfHorz_d_z_noe_nosemiz_with2A(Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(3): two endogenous states (figs 17-18)

%% ===== doPart(4): two-endo cross-tests: second state ignored =====
if doPart(4)==1
    fprintf('\n===== doPart(4): two-endo cross-tests: second state ignored =====\n')
    %% Cross-tests for two endogenous states: make the second endogenous state do nothing (the return fn
    % ignores it), so the model is the same as the one endogenous state model, and check we get the same
    % V, policy, stationary dist and stats. Done both without and with the grid interpolation layer,
    % so this checks the 2A code paths (incl. GI2A) against the one endogenous state code paths.
    output=CoreInfHorz_CrossTests_nod_nosemiz_2Aignored(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreInfHorz_CrossTests_d_nosemiz_2Aignored(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(4): two-endo cross-tests: second state ignored

%% ===== doPart(5): two-endo cross-tests: first state ignored =====
if doPart(5)==1
    fprintf('\n===== doPart(5): two-endo cross-tests: first state ignored =====\n')
    %% Same again, but now it is the FIRST endogenous state that is ignored (and the second is the actual asset)
    % Only done without the grid interpolation layer, as GI interpolates the first endogenous state,
    % which here is the ignored dummy, so a GI version would not test anything meaningful.
    output=CoreInfHorz_CrossTests_nod_nosemiz_1Aignored(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreInfHorz_CrossTests_d_nosemiz_1Aignored(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(5): two-endo cross-tests: first state ignored

%% Done!

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreInfHorzTestsdiary.txt')

diary off
