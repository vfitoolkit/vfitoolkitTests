% PType tests for the Core FHorz commands.
%
% Tests 1-4 use only z (no e, no semiz) and rely on existing return functions
% from CoreFHorzTests/CoreFHorz_ReturnFns/.
%
% The ShockTests test at the end uses all 8 (z,e,semiz) combinations across 8
% PTypes and must be set up via Names_i + per-type structures, because the
% n_z, z_grid, pi_z and vfoptions pieces differ across types.
%
% No figures are drawn anywhere in this bank, so only the diary is saved.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): N_i numeric vs Names_i cell
% doPart(2): N_i vs no-PType with identity z
% doPart(3): N_i all types identical vs one solve
% doPart(4): N_i two distinct types vs stacked solves
% doPart(5): jequaloneDist accepted three ways: jequaloneDist accepted three ways
% doPart(6): per-type N_j, different lifespans
% doPart(7): ShockTests: ShockTests: 8 shock combos
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1];

%% Diary of the command window output
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzPTypeTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzPTypeTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzPTypeTestsdiary.txt
fprintf('CoreFHorzPTypeTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreFHorzPTypeTests_subcodes/')
addpath('./CoreFHorzPTypeTests_subcodes/ShockTests/')
addpath('./CoreFHorzPTypeTests_Setup/')

% Reuse return functions and the existing setup from CoreFHorzTests
addpath('../CoreFHorzTests/CoreFHorz_ReturnFns/')
addpath('../CoreFHorzTests/CoreFHorz_ReturnFns/Semiz_ReturnFns/')
addpath('../CoreFHorzTests/CoreFHorzTests_Setup/')

%% Setup: builds on CoreFHorz_setup and adds PType-specific pieces
CoreFHorzPType_setup

%% ===== doPart(1): N_i numeric vs Names_i cell =====
if doPart(1)==1
    fprintf('\n===== doPart(1): N_i numeric vs Names_i cell =====\n')
    %% 1. N_i (numeric) vs Names_i (cell) — should give identical output
    output=CoreFHorzPType_NivsNames(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,Names_i);
end % doPart(1): N_i numeric vs Names_i cell

%% ===== doPart(2): N_i vs no-PType with identity z =====
if doPart(2)==1
    fprintf('\n===== doPart(2): N_i vs no-PType with identity z =====\n')
    %% 2. N_i vs no-PType with z having identity transitions (z encodes the type)
    output=CoreFHorzPType_NivsZidentity(n_d,n_a,N_j,d_grid,a_grid,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i);
end % doPart(2): N_i vs no-PType with identity z

%% ===== doPart(3): N_i all types identical vs one solve =====
if doPart(3)==1
    fprintf('\n===== doPart(3): N_i all types identical vs one solve =====\n')
    %% 3. N_i with all types identical ≡ a single no-PType solve
    output=CoreFHorzPType_NiAllIdentical_vsNoPType(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i);
end % doPart(3): N_i all types identical vs one solve

%% ===== doPart(4): N_i two distinct types vs stacked solves =====
if doPart(4)==1
    fprintf('\n===== doPart(4): N_i two distinct types vs stacked solves =====\n')
    %% 4. N_i with two distinct types ≡ two no-PType solves stacked
    output=CoreFHorzPType_NiTwoTypes_vsTwoSolves(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames);
end % doPart(4): N_i two distinct types vs stacked solves

%% ===== doPart(5): jequaloneDist accepted three ways: jequaloneDist accepted three ways =====
if doPart(5)==1
    fprintf('\n===== doPart(5): jequaloneDist accepted three ways: jequaloneDist accepted three ways =====\n')
    %% 5. jequaloneDist accepted three ways: [n_a,n_z], [n_a,n_z,N_i], struct
    output=CoreFHorzPType_jequaloneDist_3ways(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,Names_i);
end % doPart(5): jequaloneDist accepted three ways: jequaloneDist accepted three ways

%% ===== doPart(6): per-type N_j, different lifespans =====
if doPart(6)==1
    fprintf('\n===== doPart(6): per-type N_j, different lifespans =====\n')
    %% 6. Per-type N_j (different lifespans) — PType as struct ≡ stacked solo solves
    output=CoreFHorzPType_PerTypeNj(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames);
end % doPart(6): per-type N_j, different lifespans

%% ===== doPart(7): ShockTests: ShockTests: 8 shock combos =====
if doPart(7)==1
    fprintf('\n===== doPart(7): ShockTests: ShockTests: 8 shock combos =====\n')
    %% ShockTests: one PType solve with 8 different (z,e,semiz) combos
    % Requires Names_i (per-type structures) because n_z/z_grid/pi_z and the e/semiz
    % pieces of vfoptions differ across types.
    output=CoreFHorzPType_ShockTests_8types(n_d,n_a,n_z,n_d_semiz,d_grid_semiz,n_d2_semiz,d2_grid_semiz,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(7): ShockTests: ShockTests: 8 shock combos

% All looks good!

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreFHorzPTypeTestsdiary.txt')

diary off
