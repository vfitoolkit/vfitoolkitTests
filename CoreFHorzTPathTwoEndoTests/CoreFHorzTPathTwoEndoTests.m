% Implement core TPath tests of the VFI Toolkit FHorz commands with two endogenous states.
% with/without d
% with/without z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% with/without semiz [NOT YET IMPLEMENTED]


%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): nosemiz (figs 1-8)
%
% One part for now, because this bank is currently a single tier of 8 variants, so doPart is a
% scalar (checkcode objects to a one-element []). doPart(1) indexes a scalar exactly as it indexes
% a vector, so when the next tier is added this becomes doPart=[1,1]; and nothing else changes.
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=1;

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzTPathTwoEndoTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzTPathTwoEndoTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzTPathTwoEndoTestsdiary.txt
fprintf('CoreFHorzTPathTwoEndoTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

%%
addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreFHorzTPathTwoEndoTests_subcodes/')
addpath('./CoreFHorzTPathTwoEndoTests_Setup/')
addpath('./CoreFHorzTPathTwoEndo_ReturnFns/')

%% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorzTPathTwoEndo_setup

%% ===== doPart(1): nosemiz (figs 1-8) =====
if doPart(1)==1
    fprintf('\n===== doPart(1): nosemiz (figs 1-8) =====\n')
    %% without d, without z, without e, without semiz
    figure_c=1;
    output=CoreFHorzTPathTwoEndo_nod_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTwoEndoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, without z, without e, without semiz
    figure_c=2;
    output=CoreFHorzTPathTwoEndo_d_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTwoEndoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, with z, without e, without semiz
    figure_c=3;
    output=CoreFHorzTPathTwoEndo_nod_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTwoEndoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z, without e, without semiz
    figure_c=4;
    output=CoreFHorzTPathTwoEndo_d_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTwoEndoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, without z, with e, without semiz
    figure_c=5;
    output=CoreFHorzTPathTwoEndo_nod_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTwoEndoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, without z, with e, without semiz
    figure_c=6;
    output=CoreFHorzTPathTwoEndo_d_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTwoEndoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, with z, with e, without semiz
    figure_c=7;
    output=CoreFHorzTPathTwoEndo_nod_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTwoEndoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z, with e, without semiz
    figure_c=8;
    n_a_notsobig=[251,4];
    a_grid_notsobig=[a_grid_big(round(linspace(1,n_a_big(1),n_a_notsobig(1)))'); a_grid_big(n_a_big(1)+1:end)];
    output=CoreFHorzTPathTwoEndo_d_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_notsobig,n_z,N_j,d_grid,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTwoEndoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(1): nosemiz (figs 1-8)

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreFHorzTPathTwoEndoTestsdiary.txt')

diary off
