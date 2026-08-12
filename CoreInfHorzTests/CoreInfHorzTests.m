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
%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreInfHorzTestsdiary.txt','file')
    delete('./TestOutput/CoreInfHorzTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreInfHorzTestsdiary.txt


addpath('./CoreInfHorzTests_subcodes/')
addpath('./CoreInfHorzTests_Setup/')
addpath('./CoreInfHorz_ReturnFns/')
addpath('./CoreInfHorzTests_subcodes/CrossTests/')
% Setup so that use the same d,a,z in all the models that use them
CoreInfHorz_setup

%% without d, without z
figure_c=1;
output=CoreInfHorz_nod_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% Note: no exogenous shock, so the stationary dist is a single mass point and
% the distribution statistics are degenerate. Solver/GI/ValueFnFromPolicy tests
% are the meaningful ones here.

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

%% Cross-tests: run things that should give the same answer, disguised different ways, and check they do
% A: a single trivial z point (value 1, prob 1) reproduces the noz code path
% B: an iid disguised as a markov (identical rows) => stationary marginal over z equals that row
% C: relabelling (permuting) the z states leaves all aggregate statistics unchanged
output=CoreInfHorz_CrossTests_nod_nosemiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);

output=CoreInfHorz_CrossTests_d_nosemiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline);

%% Two endogenous states (mirror figs 3 & 4, but with a second endogenous state)
% Model: simplified Kitao (2008), no taxes -- assets + occupation (worker/entrepreneur).
% fig 12 adds Bruggemann (2021) endogenous labor supply as the decision variable d.

%% without d, with z, two endogenous states
figure_c=11;
output=CoreInfHorz_nod_z_noe_nosemiz_with2A(Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d (endogenous labor), with z, two endogenous states
figure_c=12;
output=CoreInfHorz_d_z_noe_nosemiz_with2A(Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% Done!

diary off
