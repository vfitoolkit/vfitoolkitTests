% Implement lots of tests of the VFI Toolkit FHorz transition-path commands
% under Quasi-Hyperbolic discounting (Naive and Sophisticated)
% with/without d
% with/without z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
% with/without fastOLG
%
% with/without semiz [NOT SUPPORTED by ValueFnOnTransPath_FHorz_QuasiHyperbolic]
% experienceasset    [NOT SUPPORTED by ValueFnOnTransPath_FHorz_QuasiHyperbolic]


%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreFHorzTPathQHTestsdiary.txt','file')
    delete('../TestOutput/CoreFHorzTPathQHTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreFHorzTPathQHTestsdiary.txt

%%
addpath('../CoreFHorzTPathTests_Setup/')
addpath('../CoreFHorzTPath_ReturnFns/')

addpath('./withQuasiHyperbolicDiscounting_subcodes/')
addpath('./withQuasiHyperbolicDiscounting_subcodes/CrossTests/')
% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorzTPath_setup

% Add the QH-specific bits to the baseline so each subcode picks them up via vfoptionsbaseline
Params.beta0=0.9; % additional today-tomorrow discount factor
vfoptionsbaseline.QHadditionaldiscount={'beta0'};
% vfoptions.exoticpreferences='QuasiHyperbolic';
% vfoptions.quasi_hyperbolic='Naive';        % set inside each subcode
% vfoptions.quasi_hyperbolic='Sophisticated';

% Note: for QH, V_final passed to ValueFnOnTransPath_Case1_FHorz is the 3rd output of
% ValueFnIter_Case1_FHorz (Valt for Naive, Vunderbar for Sophisticated). Subcodes build it from the
% steady-state solve. [robertdkirkby: add warnings/errors at toolkit boundary later.]

%% Uncomment to test with age-dependent shocks
% z_grid=z_grid.*ones(1,N_j,'gpuArray');
% pi_z=pi_z.*ones(1,1,N_j,'gpuArray');
% vfoptionsbaseline.pi_e=vfoptionsbaseline.pi_e.*ones(1,N_j,'gpuArray');
% vfoptionsbaseline.e_grid=vfoptionsbaseline.e_grid.*ones(1,N_j,'gpuArray');
% simoptionsbaseline.pi_e=simoptionsbaseline.pi_e.*ones(1,N_j,'gpuArray');
% simoptionsbaseline.pi_e=simoptionsbaseline.pi_e.*ones(1,N_j,'gpuArray');

%% without d, without z, without e, without semiz
figure_c=1;
output=QHDFHorzTPath_nod_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathQHTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% with d, without z, without e, without semiz
figure_c=2;
output=QHDFHorzTPath_d_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathQHTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% without d, with z, without e, without semiz
figure_c=3;
output=QHDFHorzTPath_nod_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathQHTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% with d, with z, without e, without semiz
figure_c=4;
output=QHDFHorzTPath_d_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathQHTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% without d, without z, with e, without semiz
figure_c=5;
output=QHDFHorzTPath_nod_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathQHTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% with d, without z, with e, without semiz
figure_c=6;
output=QHDFHorzTPath_d_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathQHTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% without d, with z, with e, without semiz
figure_c=7;
output=QHDFHorzTPath_nod_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathQHTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% with d, with z, with e, without semiz
figure_c=8;
output=QHDFHorzTPath_d_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathQHTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
% Note: likely redundant given the exponential-discounting cross-tests already pass, but test the
% unlikely case that QH (beta0~=1) introduces a difference somewhere.
output=QHDFHorzTPath_CrossTests_nod_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);

output=QHDFHorzTPath_CrossTests_d_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);

% all looking good :)


%% Done! Damn that was a lot of tests. Glad that is over.


diary off

%% THINGS NOT CHECKED
% semiz with QH on a transition path (not supported by ValueFnOnTransPath_FHorz_QuasiHyperbolic)
% experienceasset with QH on a transition path (errors out explicitly)
