% Implement tests of the core VFI Toolkit FHorz with ExpAssetsemiz commands
% experienceassetsemiz: aprime depends on (d2,a2,semiz), so semiz is always present
% (semiz drives the experience asset). semiz is the semi-exogenous state, so the
% semiz decision d3 is also always present.
% with/without d1
% with/without ordinary z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% Tier ordering follows CoreFHorzExpAssetTests and CoreFHorzExpAssetzTests:
%   noa1 (figs 1-8):    the experienceassetsemiz a2 is the only endogenous state
%   withA1 (figs 9-16): a1 (standard endogenous state) alongside the experienceassetsemiz
%   with2A1 (figs 17-24): two standard endogenous assets, a = [a1_1, a1_2 (binary), a2]

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): noa1 figs 1-8
% doPart(2): noa1 cross-tests
% doPart(3): withA1 figs 9-16
% doPart(4): withA1 cross-tests
% doPart(5): With TWO standard endogenous figs 17-24
% doPart(6): With TWO standard endogenous cross-tests
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzExpAssetsemizTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzExpAssetsemizTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzExpAssetsemizTestsdiary.txt
fprintf('CoreFHorzExpAssetsemizTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreFHorzExpAssetsemizTests_subcodes/')
addpath('./CoreFHorzExpAssetsemizTests_Setup/')
addpath('./CoreFHorzExpAssetsemiz_ReturnFns/')
addpath('./CoreFHorzExpAssetsemizTests_subcodes/CrossTests/')


%% Setup so that use the same d,a,semiz,z,e in all the models that use them
CoreFHorzExpAssetsemiz_setup


%% ================= noa1 (figs 1-8): the experience asset a2 is the ONLY endogenous state =================
% TEST-FIRST: these error when run, as the ExpAssetsemiz noa1 raws do not exist in the toolkit yet
% (ValueFnIter_FHorz_ExpAssetsemiz errors on N_a1==0).
addpath('./CoreFHorzExpAssetsemizTests_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzExpAssetsemiz_ReturnFns/Noa1_ReturnFns/')

%% ===== doPart(1): noa1 figs 1-8 =====
if doPart(1)==1
    fprintf('\n===== doPart(1): noa1 figs 1-8 =====\n')
    %% without d1, without z, without e, noa1
    figure_c=1;
    output=CoreFHorzExpAssetsemiz_nod1_noz_noe_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, without e, noa1
    figure_c=2;
    output=CoreFHorzExpAssetsemiz_d1_noz_noe_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, without z, with e, noa1
    figure_c=3;
    output=CoreFHorzExpAssetsemiz_nod1_noz_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, with e, noa1
    figure_c=4;
    output=CoreFHorzExpAssetsemiz_d1_noz_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, without e, noa1
    figure_c=5;
    output=CoreFHorzExpAssetsemiz_nod1_z_noe_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, without e, noa1
    figure_c=6;
    output=CoreFHorzExpAssetsemiz_d1_z_noe_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, with e, noa1
    figure_c=7;
    output=CoreFHorzExpAssetsemiz_nod1_z_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, with e, noa1
    figure_c=8;
    output=CoreFHorzExpAssetsemiz_d1_z_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(1): noa1 figs 1-8

%% ===== doPart(2): noa1 cross-tests =====
if doPart(2)==1
    fprintf('\n===== doPart(2): noa1 cross-tests =====\n')
    %% noa1 cross-tests (TEST-FIRST: error when run, same reason as above)
    % CrossTest 4: noa1 vs withA1 model where a1 (n_a1=1) is ignored (should match bit-exact)
    output=CoreFHorzExpAssetsemiz_CrossTests4_nod1_noz_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests4_d1_noz_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests4_nod1_z_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests4_d1_z_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(2): noa1 cross-tests


%% ================= withA1 (figs 9-16): a1 alongside the experienceassetsemiz =================

%% ===== doPart(3): withA1 figs 9-16 =====
if doPart(3)==1
    fprintf('\n===== doPart(3): withA1 figs 9-16 =====\n')
    %% without d1, without z, without e
    figure_c=9;
    output=CoreFHorzExpAssetsemiz_nod1_noz_noe(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, without e
    n_a_notsobig=[301,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=10;
    output=CoreFHorzExpAssetsemiz_d1_noz_noe(n_d_withd1,n_a,n_a_notsobig,0,N_j,d_grid_withd1,a_grid,a_grid_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, without z, with e
    n_a_notsobig=[301,13];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=11;
    output=CoreFHorzExpAssetsemiz_nod1_noz_e(n_d_withoutd1,n_a,n_a_notsobig,0,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, with e
    n_a_notsobig=[201,13];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=12;
    output=CoreFHorzExpAssetsemiz_d1_noz_e(n_d_withd1,n_a,n_a_notsobig,0,N_j,d_grid_withd1,a_grid,a_grid_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, without e
    n_a_notsobig=[301,13];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=13;
    output=CoreFHorzExpAssetsemiz_nod1_z_noe(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, without e
    n_a_notsobig=[301,13];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=14;
    output=CoreFHorzExpAssetsemiz_d1_z_noe(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, with e
    n_a_notsobig=[301,13];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=15;
    output=CoreFHorzExpAssetsemiz_nod1_z_e(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, with e
    n_a_notsobig=[201,13];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=16;
    output=CoreFHorzExpAssetsemiz_d1_z_e(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(3): withA1 figs 9-16


%% ===== doPart(4): withA1 cross-tests =====
if doPart(4)==1
    fprintf('\n===== doPart(4): withA1 cross-tests =====\n')
    %% withA1 cross-tests
    % CrossTest 1: experienceassetsemiz with a degenerate (d3-invariant) semiz vs experienceassetz (should match)
    output=CoreFHorzExpAssetsemiz_CrossTests_nod1(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests_d1(n_d_withd1,n_a,n_a_big,0,N_j,d_grid_withd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    % CrossTest 2: 'fake' experienceassetsemiz that ignores semiz vs plain experienceasset+semiz (should match)
    output=CoreFHorzExpAssetsemiz_CrossTests2_nod1(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests2_d1(n_d_withd1,n_a,n_a_big,0,N_j,d_grid_withd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    % CrossTest 3: experienceassetsemiz+z vs experienceassetz with combined z=[semiz,z] (pins bothz ordering)
    output=CoreFHorzExpAssetsemiz_CrossTests3_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests3_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(4): withA1 cross-tests


%% ================= With TWO standard endogenous assets (figs 17-24) =================
% a = [a1_1, a1_2 (binary), a2 (experienceassetsemiz)]
% TEST-FIRST: these error when run, as the ExpAssetsemiz DC2A/GI2A/DC2A_GI2A raws do not
% exist in the toolkit yet.
addpath('./CoreFHorzExpAssetsemizTests_subcodes/With2A1_subcodes/')
addpath('./CoreFHorzExpAssetsemiz_ReturnFns/With2A1_ReturnFns/')

% a1main is kept modest here because the binary second asset doubles the a-grid.
% n_a_2A1=[a1, a1_2, a2] and a_grid_2A1 come from CoreFHorzExpAssetsemiz_setup; a1_2 is a genuine
% multi-point second standard asset, so the subcodes take n_a/a_grid as given and build nothing.
n_a_2A1_notsobig=[151,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

%% ===== doPart(5): With TWO standard endogenous figs 17-24 =====
if doPart(5)==1
    fprintf('\n===== doPart(5): With TWO standard endogenous figs 17-24 =====\n')
    %% with2A1, without d1, without z, without e
    figure_c=17;
    output=CoreFHorzExpAssetsemiz_nod1_noz_noe_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A1, with d1, without z, without e
    figure_c=18;
    output=CoreFHorzExpAssetsemiz_d1_noz_noe_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A1, without d1, without z, with e
    figure_c=19;
    output=CoreFHorzExpAssetsemiz_nod1_noz_e_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A1, with d1, without z, with e
    figure_c=20;
    output=CoreFHorzExpAssetsemiz_d1_noz_e_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A1, without d1, with z, without e
    figure_c=21;
    output=CoreFHorzExpAssetsemiz_nod1_z_noe_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A1, with d1, with z, without e
    figure_c=22;
    output=CoreFHorzExpAssetsemiz_d1_z_noe_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A1, without d1, with z, with e
    figure_c=23;
    output=CoreFHorzExpAssetsemiz_nod1_z_e_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A1, with d1, with z, with e
    figure_c=24;
    output=CoreFHorzExpAssetsemiz_d1_z_e_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(5): With TWO standard endogenous figs 17-24

%% ===== doPart(6): With TWO standard endogenous cross-tests =====
if doPart(6)==1
    fprintf('\n===== doPart(6): With TWO standard endogenous cross-tests =====\n')
    %% with2A1 cross-tests (TEST-FIRST: error when run, same reason as above)
    % CrossTest 5: a degenerate second standard asset a1_2 (single point {0}) reduces the
    % two-standard-asset model back to the with-a1 model (i.e. DC2A/GI2A vs the ordinary solvers)
    output=CoreFHorzExpAssetsemiz_CrossTests5_nod1_with2A1(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests5_d1_with2A1(n_d_withd1,n_a,n_a_big,0,N_j,d_grid_withd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(6): With TWO standard endogenous cross-tests

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreFHorzExpAssetsemizTestsdiary.txt')

diary off
