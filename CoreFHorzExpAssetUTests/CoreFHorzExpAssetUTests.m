% Implement lots of tests of the core VFI Toolkit FHorz with ExpAssetU commands
% with/without d1
% with/without z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% with/without semiz
%
% This is all done WITHOUT a1 first (experience asset is the only endogenous state), then rerun
% all 16 WITH a1 (where divide-and-conquer and grid interpolation layer also apply), then again
% WITH a genuine second standard asset a1_2 (with2A1, which triggers DC2A/GI2A/DC2A_GI2A).

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1):  noa1 nosemiz (figs 1-8)
% doPart(2):  noa1 nosemiz cross-tests
% doPart(3):  noa1 semiz (figs 9-16)
% doPart(4):  noa1 semiz cross-tests
% doPart(5):  a1 nosemiz (figs 17-24)
% doPart(6):  a1 nosemiz cross-tests
% doPart(7):  a1 nosemiz fake-experienceassetu cross-tests
% doPart(8):  a1 semiz (figs 25-32)
% doPart(9):  a1 semiz cross-tests
% doPart(10): a1 semiz cross-tests 2 (semi-exo that is really a markov)
% doPart(11): a1 semiz fake-experienceassetu cross-tests
% doPart(12): 2a1 nosemiz (figs 33-40)
% doPart(13): 2a1 nosemiz cross-tests
% doPart(14): 2a1 semiz (figs 41-48) - WILL ERROR: there is no ExpAssetUSemiExo DC2A/GI2A family
%             in the toolkit yet, so this part is here test-first. It is left at 1 because that is
%             what the bank does today; set it to 0 to get a clean run of everything else.
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so Fig 42
% is the same test whatever doPart says, and a png from a previous run is never overwritten by a
% different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
%
% NOTE: figs 40, 44, 46 and 48 stay commented out INSIDE their parts - they out-of-memory on this
% machine. doPart is for whole tiers; an individual call that cannot run stays commented.
% NOTE: figs 44, 46, 47 and 48 stay commented out INSIDE part 14 - they out-of-memory on this
% machine (fig 47 crashed MATLAB outright on the 2026-10-12 run: it is z+e+semiz+2A1 on the big
% grid with DC2A+GI2A, the heaviest call in the bank). doPart is for whole tiers; an individual
% call that cannot run stays commented.
%
% This run: parts 1-4, the noa1 tiers. Their ReturnFns gained a 'transfer' parameter on
% 2026-10-12 so that working-age consumption is positive at a=0 (a is the unconsumable
% experience asset, so earnings there were exactly 0, F=-Inf, and the experience-asset
% interpolation propagated that -Inf to every state and backwards through every age, leaving
% whole blocks of this tier's checks vacuous - the 'V_Jplus1 with age-dependent shocks ... NaN'
% lines were the visible symptom). Parts 1-4 are exactly the parts that call the changed
% subcodes, so they are the ones whose numbers move; parts 5-14 are untouched and were clean on
% the 2026-10-12 runs. Parts are independent (no part reads another's output and nothing is
% passed through a .mat), so running this subset is safe.
doPart=[1,1,1,1,0,0,0,0,0,0,0,0,0,0];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzExpAssetUTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzExpAssetUTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzExpAssetUTestsdiary.txt
fprintf('CoreFHorzExpAssetUTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreFHorzExpAssetUTests_subcodes/')
addpath('./CoreFHorzExpAssetUTests_Setup/')
addpath('./CoreFHorzExpAssetU_ReturnFns/')
addpath('./CoreFHorzExpAssetUTests_subcodes/CrossTests/')
% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorzExpAssetU_setup

%% ================= WITHOUT a1 (figs 1-16): experience asset is the only endogenous state =================
% No DC/GI/DC+GI blocks (irrelevant without a1).
% Pass n_a_justexpasset as n_a, a_grid_justexpasset as a_grid. n_a_big/a_grid_big slots unused.

addpath('./CoreFHorzExpAssetUTests_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzExpAssetUTests_subcodes/Noa1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAssetU_ReturnFns/Noa1_ReturnFns/')
addpath('./CoreFHorzExpAssetU_ReturnFns/Noa1_ReturnFns/Semiz_ReturnFns/')
% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorzExpAssetU_setup

%% ===== doPart(1): noa1 nosemiz (figs 1-8) =====
if doPart(1)==1
    fprintf('\n===== doPart(1): noa1 nosemiz (figs 1-8) =====\n')
    %% noa1 nosemiz (8 variants)

    %% without d1, without z, without e, noa1, nosemiz
    figure_c=1;
    output=CoreFHorzExpAssetU_nod1_noz_noe_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, without e, noa1, nosemiz
    figure_c=2;
    output=CoreFHorzExpAssetU_d1_noz_noe_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, without e, noa1, nosemiz
    figure_c=3;
    output=CoreFHorzExpAssetU_nod1_z_noe_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, without e, noa1, nosemiz
    figure_c=4;
    output=CoreFHorzExpAssetU_d1_z_noe_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, without z, with e, noa1, nosemiz
    figure_c=5;
    output=CoreFHorzExpAssetU_nod1_noz_e_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, with e, noa1, nosemiz
    figure_c=6;
    output=CoreFHorzExpAssetU_d1_noz_e_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, with e, noa1, nosemiz
    figure_c=7;
    output=CoreFHorzExpAssetU_nod1_z_e_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, with e, noa1, nosemiz
    figure_c=8;
    output=CoreFHorzExpAssetU_d1_z_e_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(1): noa1 nosemiz (figs 1-8)

%% ===== doPart(2): noa1 nosemiz cross-tests =====
if doPart(2)==1
    fprintf('\n===== doPart(2): noa1 nosemiz cross-tests =====\n')
    %% noa1 nosemiz cross-tests
    % Markov-as-iid equivalence cross-test
    output=CoreFHorzExpAssetU_CrossTests_nod1_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetU_CrossTests_d1_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % 'fake' experienceassetu vs experienceasset
    output=CoreFHorzExpAssetU_CrossTests3_nod1_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetU_CrossTests3_d1_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % no a1 vs model with a1 but where it is ignored
    output=CoreFHorzExpAssetU_CrossTests4_nod1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetU_CrossTests4_d1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(2): noa1 nosemiz cross-tests




%% ===== doPart(3): noa1 semiz (figs 9-16) =====
if doPart(3)==1
    fprintf('\n===== doPart(3): noa1 semiz (figs 9-16) =====\n')
    %% noa1 semiz (8 variants)

    %% without d1, without z, without e, noa1, semiz
    figure_c=9;
    output=CoreFHorzExpAssetU_nod1_noz_noe_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, without e, noa1, semiz
    figure_c=10;
    output=CoreFHorzExpAssetU_d1_noz_noe_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, without e, noa1, semiz
    figure_c=11;
    output=CoreFHorzExpAssetU_nod1_z_noe_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, without e, noa1, semiz
    figure_c=12;
    output=CoreFHorzExpAssetU_d1_z_noe_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, without z, with e, noa1, semiz
    figure_c=13;
    output=CoreFHorzExpAssetU_nod1_noz_e_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, with e, noa1, semiz
    figure_c=14;
    output=CoreFHorzExpAssetU_d1_noz_e_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, with e, noa1, semiz
    figure_c=15;
    output=CoreFHorzExpAssetU_nod1_z_e_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, with e, noa1, semiz
    figure_c=16;
    output=CoreFHorzExpAssetU_d1_z_e_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(3): noa1 semiz (figs 9-16)

%% ===== doPart(4): noa1 semiz cross-tests =====
if doPart(4)==1
    fprintf('\n===== doPart(4): noa1 semiz cross-tests =====\n')
    %% noa1 semiz cross-tests
    % Markov-as-iid equivalence cross-test
    output=CoreFHorzExpAssetU_CrossTests_nod1_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetU_CrossTests_d1_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % semiz state as a plain Markov should give same answer
    output=CoreFHorzExpAssetU_CrossTests2_nod1_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetU_CrossTests2_d1_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % 'fake' experienceassetu vs experienceasset
    output=CoreFHorzExpAssetU_CrossTests3_nod1_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetU_CrossTests3_d1_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % no a1 vs model with a1 but where it is ignored
    output=CoreFHorzExpAssetU_CrossTests4_nod1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetU_CrossTests4_d1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(4): noa1 semiz cross-tests







%% ================= WITH a1 (figs 17-32) =================
% Reset the setup so this half can be run independently of the WITHOUT-a1 half above.
CoreFHorzExpAssetU_setup

%% ===== doPart(5): a1 nosemiz (figs 17-24) =====
if doPart(5)==1
    fprintf('\n===== doPart(5): a1 nosemiz (figs 17-24) =====\n')
    %% without d1, without z, without e, without semiz
    figure_c=17;
    output=CoreFHorzExpAssetU_nod1_noz_noe_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d1, without z, without e, without semiz
    figure_c=18;
    n_a_notsobig=[501,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    output=CoreFHorzExpAssetU_d1_noz_noe_nosemiz(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% without d1, with z, without e, without semiz
    figure_c=19;
    n_a_notsobig=[501,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    output=CoreFHorzExpAssetU_nod1_z_noe_nosemiz(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% with d1, with z, without e, without semiz
    n_a_notsobig=[201,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=20;
    output=CoreFHorzExpAssetU_d1_z_noe_nosemiz(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% without d1, without z, with e, without semiz
    figure_c=21;
    n_a_notsobig=[501,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    output=CoreFHorzExpAssetU_nod1_noz_e_nosemiz(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% with d1, without z, with e, without semiz
    n_a_notsobig=[201,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=22;
    output=CoreFHorzExpAssetU_d1_noz_e_nosemiz(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% without d1, with z, with e, without semiz
    n_a_notsobig=[301,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=23;
    output=CoreFHorzExpAssetU_nod1_z_e_nosemiz(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% with d1, with z, with e, without semiz
    n_a_notsobig=[201,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=24;
    output=CoreFHorzExpAssetU_d1_z_e_nosemiz(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good (as good as it can be expected to given the n_a_notsobig)
end % doPart(5): a1 nosemiz (figs 17-24)

%% ===== doPart(6): a1 nosemiz cross-tests =====
if doPart(6)==1
    fprintf('\n===== doPart(6): a1 nosemiz cross-tests =====\n')
    %% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
    output=CoreFHorzExpAssetU_CrossTests_nod1_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorzExpAssetU_CrossTests_d1_nosemiz(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % all looking good :)
end % doPart(6): a1 nosemiz cross-tests

%% ===== doPart(7): a1 nosemiz fake-experienceassetu cross-tests =====
if doPart(7)==1
    fprintf('\n===== doPart(7): a1 nosemiz fake-experienceassetu cross-tests =====\n')
    %% Do a test with a 'fake experienceassetu' and compare to experienceasset
    output=CoreFHorzExpAssetU_CrossTests3_nod1_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorzExpAssetU_CrossTests3_d1_nosemiz(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % all looking good :)
end % doPart(7): a1 nosemiz fake-experienceassetu cross-tests









%% Worth doing a 'clear all' here, but not necessary.
% Mainly is so you can run second half independent of first half

%% That is all the without semiz, now with semiz
% From here on, it is the eight with semiz
% From here on, use n_d_semiz and d_grid_semiz as the inputs (instead of n_d and d_grid)

addpath('./CoreFHorzExpAssetUTests_subcodes/')
addpath('./CoreFHorzExpAssetUTests_Setup/')
addpath('./CoreFHorzExpAssetU_ReturnFns/')
addpath('./CoreFHorzExpAssetUTests_subcodes/CrossTests/')

addpath('./CoreFHorzExpAssetUTests_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAssetU_ReturnFns/Semiz_ReturnFns/')
% Uses the same setup, which already had a semi-exogenous state, just that it wasn't used.
CoreFHorzExpAssetU_setup

% For models without d1, use:
% n_d2_semiz and d2_grid_semiz (as n_d and d_grid)
% For models with d1, use:
% n_d_semiz and d_grid_semiz (as n_d and d_grid)

%% ===== doPart(8): a1 semiz (figs 25-32) =====
if doPart(8)==1
    fprintf('\n===== doPart(8): a1 semiz (figs 25-32) =====\n')
    %% without d1, without z, without e, with semiz
    figure_c=25;
    n_a_notsobig=[201,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    output=CoreFHorzExpAssetU_nod1_noz_noe_semiz(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% with d1, without z, without e, with semiz
    n_a_notsobig=[201,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=26;
    output=CoreFHorzExpAssetU_d1_noz_noe_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% without d1, with z, without e, with semiz
    n_a_notsobig=[301,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=27;
    output=CoreFHorzExpAssetU_nod1_z_noe_semiz(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% with d1, with z, without e, with semiz
    n_a_notsobig=[151,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=28;
    output=CoreFHorzExpAssetU_d1_z_noe_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% without d1, without z, with e, with semiz
    n_a_notsobig=[301,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=29;
    output=CoreFHorzExpAssetU_nod1_noz_e_semiz(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% with d1, without z, with e, with semiz
    n_a_notsobig=[201,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=30;
    output=CoreFHorzExpAssetU_d1_noz_e_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% without d1, with z, with e, with semiz
    n_a_notsobig=[151,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=31;
    output=CoreFHorzExpAssetU_nod1_z_e_semiz(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% with d1, with z, with e, with semiz
    n_a_notsobig=[151,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=32;
    output=CoreFHorzExpAssetU_d1_z_e_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    % looks good (as good as it can be expected to given the n_a_notsobig)
end % doPart(8): a1 semiz (figs 25-32)

%% ===== doPart(9): a1 semiz cross-tests =====
if doPart(9)==1
    fprintf('\n===== doPart(9): a1 semiz cross-tests =====\n')
    %% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
    output=CoreFHorzExpAssetU_CrossTests_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorzExpAssetU_CrossTests_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % all looking good :)
end % doPart(9): a1 semiz cross-tests

%% ===== doPart(10): a1 semiz cross-tests 2 (semi-exo that is really a markov) =====
if doPart(10)==1
    fprintf('\n===== doPart(10): a1 semiz cross-tests 2 (semi-exo that is really a markov) =====\n')
    %% Now some further cross-tests, using a semi-exo that is really just a markov
    output=CoreFHorzExpAssetU_CrossTests2_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorzExpAssetU_CrossTests2_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % all looking good :)
end % doPart(10): a1 semiz cross-tests 2 (semi-exo that is really a markov)

%% ===== doPart(11): a1 semiz fake-experienceassetu cross-tests =====
if doPart(11)==1
    fprintf('\n===== doPart(11): a1 semiz fake-experienceassetu cross-tests =====\n')
    %% Do a test with a 'fake experienceassetu' and compare to experienceasset
    output=CoreFHorzExpAssetU_CrossTests3_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorzExpAssetU_CrossTests3_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % all looking good :)
end % doPart(11): a1 semiz fake-experienceassetu cross-tests




















%% ================= WITH 2a1 (figs 33-48): two standard endogenous assets (triggers DC2A/GI2A/DC2A_GI2A) =================
% A genuine multi-point second standard asset a1_2 (return r2) is spliced between the liquid asset
% a1 and the experience asset a2: n_a_2A1=[n_a1, n_a1_2, n_a2] (built in CoreFHorzExpAssetU_setup).
% Two standard assets -> length(n_a1)>1 -> DC2A / GI2A / DC2A_GI2A.
% (ExpAssetU also carries the u dimension, so the notsobig a1 grids are kept conservative.)
% NOTE: the 8 semiz variants (figs 41-48) will ERROR when run -- there is no ExpAssetUSemiExo
% DC2A/GI2A family in the toolkit yet; they are included test-first for when that family is built.
CoreFHorzExpAssetU_setup
addpath('./CoreFHorzExpAssetUTests_subcodes/With2A1_subcodes/')
addpath('./CoreFHorzExpAssetUTests_subcodes/With2A1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAssetU_ReturnFns/With2A1_ReturnFns/')
addpath('./CoreFHorzExpAssetU_ReturnFns/With2A1_ReturnFns/Semiz_ReturnFns/')

%% ===== doPart(12): 2a1 nosemiz (figs 33-40) =====
if doPart(12)==1
    fprintf('\n===== doPart(12): 2a1 nosemiz (figs 33-40) =====\n')
    %% 2a1 nosemiz (8 variants)
    n_a_2A1_notsobig=[151,n_a1_2,n_a_justexpasset];
    a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
    a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

    figure_c=33;
    output=CoreFHorzExpAssetU_nod1_noz_noe_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=34;
    output=CoreFHorzExpAssetU_d1_noz_noe_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    n_a_2A1_notsobig=[101,n_a1_2,n_a_justexpasset];
    a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
    a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

    figure_c=35;
    output=CoreFHorzExpAssetU_nod1_z_noe_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=36;
    output=CoreFHorzExpAssetU_d1_z_noe_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=37;
    output=CoreFHorzExpAssetU_nod1_noz_e_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=38;
    output=CoreFHorzExpAssetU_d1_noz_e_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    n_a_2A1_notsobig=[75,n_a1_2,n_a_justexpasset];
    a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
    a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

    figure_c=39;
    output=CoreFHorzExpAssetU_nod1_z_e_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=40;
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    % COMMENTED OUT: out-of-memories on this machine, and MATLAB halts the whole script at
    % the first error, which would stop the run before it ever reaches figs 41-48.
    % output=CoreFHorzExpAssetU_d1_z_e_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
end % doPart(12): 2a1 nosemiz (figs 33-40)

%% ===== doPart(13): 2a1 nosemiz cross-tests =====
if doPart(13)==1
    fprintf('\n===== doPart(13): 2a1 nosemiz cross-tests =====\n')
    %% 2a1 nosemiz cross-tests: a degenerate second asset a1_2 (single point {0}) reduces to the with-a1 model
    output=CoreFHorzExpAssetU_CrossTests_nod1_nosemiz_with2A1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetU_CrossTests_d1_nosemiz_with2A1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(13): 2a1 nosemiz cross-tests

%% ===== doPart(14): 2a1 semiz (figs 41-48) =====
if doPart(14)==1
    fprintf('\n===== doPart(14): 2a1 semiz (figs 41-48) =====\n')
    %% 2a1 semiz (8 variants)
    n_a_2A1_notsobig=[151,n_a1_2,n_a_justexpasset];
    a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
    a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

    figure_c=41;
    output=CoreFHorzExpAssetU_nod1_noz_noe_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=42;
    output=CoreFHorzExpAssetU_d1_noz_noe_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    n_a_2A1_notsobig=[101,n_a1_2,n_a_justexpasset];
    a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
    a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

    figure_c=43;
    output=CoreFHorzExpAssetU_nod1_z_noe_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    n_a_2A1_notsobig=[75,n_a1_2,n_a_justexpasset];
    a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
    a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

    figure_c=44;
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    % output=CoreFHorzExpAssetU_d1_z_noe_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    figure_c=45;
    output=CoreFHorzExpAssetU_nod1_noz_e_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=46;
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    % output=CoreFHorzExpAssetU_d1_noz_e_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    figure_c=47;
    % output=CoreFHorzExpAssetU_nod1_z_e_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    figure_c=48;
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    % output=CoreFHorzExpAssetU_d1_z_e_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetUTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
end % doPart(14): 2a1 semiz (figs 41-48)


%% One verdict for the whole run
% The bank prints over a thousand checks; this reads the diary back and says plainly whether the
% run passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreFHorzExpAssetUTestsdiary.txt')

diary off

%% THINGS NOT CHECKED
% Check using two decision variables in any of d1 or d3 (the decision variables that are not in experience asset)
