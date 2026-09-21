% Implement core TPath tests of the VFI Toolkit FHorz commands with ExperienceAssetz.
% experienceassetz: aprime depends on (d2,a2,z), so z is always present
% with/without d1
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% with/without semiz [NOT YET IMPLEMENTED]
% noa1 [NOT YET IMPLEMENTED]
%
% This is all done with a1 (standard endogenous state alongside the experienceassetz)
% experienceassetz without a1 is not supported by VFI Toolkit


%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): figs 1-4
% doPart(2): cross-tests
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzTPathExpAssetzTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzTPathExpAssetzTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzTPathExpAssetzTestsdiary.txt
fprintf('CoreFHorzTPathExpAssetzTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

%%
addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreFHorzTPathExpAssetzTests_subcodes/')
addpath('./CoreFHorzTPathExpAssetzTests_Setup/')
addpath('./CoreFHorzTPathExpAssetz_ReturnFns/')
addpath('./CoreFHorzTPathExpAssetzTests_subcodes/CrossTests/')
% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorzTPathExpAssetz_setup

%% ===== doPart(1): figs 1-4 =====
if doPart(1)==1
    fprintf('\n===== doPart(1): figs 1-4 =====\n')
    %% without d1, with z, without e, without semiz
    figure_c=1;
    output=CoreFHorzTPathExpAssetz_nod1_z_noe_nosemiz(T,PricePath,ParamPath,n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % look good

    %% with d1, with z, without e, without semiz
    figure_c=2;
    % Following the TPath ExpAsset tests, use a smaller 'big' grid to avoid out-of-memory
    n_a_notsobig=[501,n_a_justexpasset]; % to test Grid Interpolation
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];
    output=CoreFHorzTPathExpAssetz_d1_z_noe_nosemiz(T,PricePath,ParamPath,n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % Some of the lowmemory are not quite right, seems to be just the Policy (V
    % seems fine); at first I thought L2flag, but appears to impact the DC
    % (without GI) so that is not the reason.

    %% without d1, with z, with e, without semiz
    figure_c=3;
    % experienceassetz machinery carries z-dependence, so use the smaller 'big' grid to avoid out-of-memory
    n_a_notsobig=[301,n_a_justexpasset]; % to test Grid Interpolation
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];
    output=CoreFHorzTPathExpAssetz_nod1_z_e_nosemiz(T,PricePath,ParamPath,n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % look good

    %% with d1, with z, with e, without semiz
    figure_c=4;
    n_a_notsobig=[301,n_a_justexpasset]; % to test Grid Interpolation
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];
    output=CoreFHorzTPathExpAssetz_d1_z_e_nosemiz(T,PricePath,ParamPath,n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % Some of the lowmemory are not quite right, seems to be just the Policy (V
    % seems fine); at first I thought L2flag, but appears to impact the DC
    % (without GI) so that is not the reason.
end % doPart(1): figs 1-4

%% ===== doPart(2): cross-tests =====
if doPart(2)==1
    fprintf('\n===== doPart(2): cross-tests =====\n')
    %% Cross-tests
    % CrossTest 2: 'fake' experienceassetz that ignores z vs plain experienceasset (should match)
    output=CoreFHorzTPathExpAssetz_CrossTest2_nod1(T,PricePath,ParamPath,n_d_withoutd1,n_a,n_z,N_j,d_grid_withoutd1,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline);
    output=CoreFHorzTPathExpAssetz_CrossTest2_d1(T,PricePath,ParamPath,n_d_withd1,n_a,n_z,N_j,d_grid_withd1,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline);
end % doPart(2): cross-tests

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreFHorzTPathExpAssetzTestsdiary.txt')

diary off
