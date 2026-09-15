function output=CoreInfHorzTPathAlgo_crosstest_jointdgrid(T,PricePath,ParamPath,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames)
% Cross test: a stacked-column d_grid against a joint grid with the never-chosen rows removed.
%
% The return function gives F=-Inf whenever d1+d2>1, so those combinations of the two decision
% variables can never be optimal. Dropping them from the grid must therefore leave the answer exactly
% unchanged, not merely close, which is what makes this a sharp test rather than an approximate one.
%
% The reduced grid is passed as n_d=[M,1] with an M-by-2 joint grid, where M is the number of
% combinations kept. TransitionPath_InfHorz accepts a joint grid of size prod(n_d)-by-length(n_d),
% and [prod([M,1]),length([M,1])]=[M,2], so this needs no change to the toolkit. If that turns out to
% mislead anything downstream (the policy is un-kronned using n_d), the exact-zero comparison below is
% what will catch it.

%% Two decision variables, five points each
n_d1=5; n_d2=5;
n_d=[n_d1,n_d2];
d1_grid=linspace(0,1,n_d1)';
d2_grid=linspace(0,1,n_d2)';
d_grid=[d1_grid; d2_grid]; % stacked-column form

ReturnFn=@(d1,d2,aprime,a,z,r,w,sigma) CoreInfHorzTPathAlgo_jointdgrid_ReturnFn(d1,d2,aprime,a,z,r,w,sigma);

% Built from scratch rather than from the baselines, which carry n_e=3. This model has no e variable,
% the same reason the noe subcodes in the main bank do this.
vfoptions1=struct();
simoptions1=struct();

%% A general eqm eqn that actually responds to the decision variables, so the one-iteration price
% update depends on the policy. A price update that ignored the policy would compare equal whatever
% the grids did.
FnsToEvaluate.K=@(d1,d2,aprime,a,z) a;
GeneralEqmEqns.capital=@(r,K) K;

transpathoptions=struct();
transpathoptions.maxiter=1; % one iteration is enough: the two grids either give the same policy or they do not
transpathoptions.GEnewprice=3;
transpathoptions.GEnewprice3.howtoupdate={'capital','r',0,0.1};

%% Solve it with the stacked-column d_grid
[V_final,Policy_final]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
AgentDist_initial=StationaryDist_InfHorz(Policy_final,n_d,n_a,n_z,pi_z,simoptions1,Params,[]);
PricePathStacked=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions1, vfoptions1, []);

%% Build the joint grid and drop the rows the return function makes infeasible
d_gridvals=CreateGridvals(n_d,d_grid,1); % [prod(n_d), 2], every combination
keep=(d_gridvals(:,1)+d_gridvals(:,2)<=1); % exactly the rows the return function does not send to -Inf
d_gridvals_reduced=d_gridvals(keep,:);
M=size(d_gridvals_reduced,1);
n_d_reduced=[M,1];
fprintf('Joint d_grid cross test: %i of %i combinations kept (the rest are -Inf in the return function) \n',M,size(d_gridvals,1))

%% Solve the same problem with the reduced joint grid
[V_final_j,Policy_final_j]=ValueFnIter_InfHorz(n_d_reduced,n_a,n_z,d_gridvals_reduced,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
AgentDist_initial_j=StationaryDist_InfHorz(Policy_final_j,n_d_reduced,n_a,n_z,pi_z,simoptions1,Params,[]);
PricePathJoint=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final_j, AgentDist_initial_j, n_d_reduced, n_a, n_z, d_gridvals_reduced,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions1, vfoptions1, []);

%% The dropped rows could never have been chosen, so the two must agree exactly
fprintf('Stacked d_grid vs reduced joint d_grid, V, this should be zero: %.3e \n',max(abs(V_final(:)-V_final_j(:))))
fprintf('Stacked d_grid vs reduced joint d_grid, one iteration of the TPath, this should be zero, r: %.3e \n',max(abs(PricePathJoint.r-PricePathStacked.r)))

output=1;

end
