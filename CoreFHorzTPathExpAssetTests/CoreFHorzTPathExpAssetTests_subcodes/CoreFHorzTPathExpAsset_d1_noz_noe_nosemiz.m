function output=CoreFHorzTPathExpAsset_d1_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c)

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
% ExperienceAsset (a2)
vfoptions.experienceasset=vfoptionsbaseline.experienceasset;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.experienceasset=simoptionsbaseline.experienceasset;
simoptions.aprimeFn=simoptionsbaseline.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;

n_z=0; z_grid=[]; pi_z=[];

% zero assets, mid points for any shocks
jequaloneDist_big=zeros([n_a_big],'gpuArray');
jequaloneDist_big(1,1)=1;
jequaloneDist=zeros([n_a],'gpuArray');
jequaloneDist(1,1)=1;

ReturnFn=@(d1,d2,a1prime,a1,a2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_noz_noe_nosemiz(d1,d2,a1prime,a1,a2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(d1,d2,a1prime,a1,a2) a1+a2;
FnsToEvaluate.earnings=@(d1,d2,a1prime,a1,a2,w,kappa_j) w*kappa_j*d1*d2*a2;


% Need period T for V and Policy
V_final=zeros([n_a,N_j],'gpuArray');
Policy_final=ones([3,n_a,N_j],'gpuArray');
Policy_final_GI=ones([5,n_a,N_j],'gpuArray');
% big versions of them
V_final_big=zeros([n_a_big,N_j],'gpuArray');
Policy_final_big=ones([3,n_a_big,N_j],'gpuArray');
Policy_final_big_GI=ones([5,n_a_big,N_j],'gpuArray');

%% With vs Without fastOLG
transpathoptionsslow.fastOLG=0;
vfoptions1=vfoptions;
simoptions1=simoptions;
[VPath1slow,PolicyPath1slow]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions1);

[VPath1fast,PolicyPath1fast]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);

fprintf('fastOLG, this should be zero: %2.8f \n',max(abs(VPath1slow(:)-VPath1fast(:))))
fprintf('fastOLG, this should be zero: %2.8f \n',max(abs(PolicyPath1slow(:)-PolicyPath1fast(:))))

clear VPath1fast VPath1slow PolicyPath1fast PolicyPath1slow

%% With and without divide-and-conquer (both with slowOLG)
[VPath1,PolicyPath1]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions1);

% PolicyVals1=PolicyInd2Val_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions1);

% Solve with divide-and-conquer, should give same answer
vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
simoptions2=simoptions;
[VPath2,PolicyPath2]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions2);

fprintf('Divide-and-conquer (slowOLG), this should be zero: %2.8f \n',max(abs(VPath1(:)-VPath2(:))))
fprintf('Divide-and-conquer (slowOLG), this should be zero: %2.8f \n',max(abs(PolicyPath1(:)-PolicyPath2(:))))

% Trying to understand the difference of 2 in the PolicyPath1-PolicyPath2
% Claude asked for all the following, and based on that concluded it was
% about how with and without DC meant that Policy is different for
% decisions where the agent is indifferent. I am not yet 100% convinced, so
% leaving this here until I am.
m = squeeze(any(abs(PolicyPath1 - PolicyPath2) > 0, 1));   % [N_a, N_j, T]
[a_idx, j_idx, t_idx] = ind2sub(size(m), find(m, 1, 'first'));
fprintf('First diff at (a=%d, j=%d, t=%d):\n  nonDC=[%s]  DC=[%s]\n', ...
    a_idx, j_idx, t_idx, ...
    num2str(PolicyPath1(:,a_idx,j_idx,t_idx)'), ...
    num2str(PolicyPath2(:,a_idx,j_idx,t_idx)'));
fprintf('  V: nonDC=%.15g  DC=%.15g  (eq=%d)\n', ...
    VPath1(a_idx,j_idx,t_idx), VPath2(a_idx,j_idx,t_idx), ...
    VPath1(a_idx,j_idx,t_idx)==VPath2(a_idx,j_idx,t_idx));

fprintf('\nDiff counts by j (age):\n');
for jj=1:size(m,2)
    c = nnz(m(:,jj,:));
    if c>0, fprintf('  j=%2d: %d\n', jj, c); end
end
fprintf('\nDiff counts by t:\n');
for tt=1:size(m,3)
    c = nnz(m(:,:,tt));
    if c>0, fprintf('  t=%3d: %d\n', tt, c); end
end

% Distinct (d2_nonDC, d2_DC) value pairs at diffs:
diffmask = squeeze(PolicyPath1(2,:,:,:)) ~= squeeze(PolicyPath2(2,:,:,:));
pairs = [reshape(squeeze(PolicyPath1(2,:,:,:)),[],1) reshape(squeeze(PolicyPath2(2,:,:,:)),[],1)];
pairs = pairs(diffmask(:),:);
[uniq, ~, ic] = unique(pairs, 'rows');
counts = accumarray(ic, 1);
fprintf('\n(d2_nonDC, d2_DC) -> count:\n');
for k=1:size(uniq,1)
    fprintf('  (%d, %d) -> %d\n', uniq(k,1), uniq(k,2), counts(k));
end


m = squeeze(any(abs(PolicyPath1 - PolicyPath2) > 0, 1));   % [N_a, N_j, T]

% How many j have diffs, and which?
fprintf('j values with diffs: ');
for jj=1:size(m,2)
    if nnz(m(:,jj,:)) > 0, fprintf('%d ', jj); end
end
fprintf('\n');

% How many t have diffs?
fprintf('t values with diffs (count): %d distinct\n', ...
        nnz(squeeze(sum(sum(m,1),2)) > 0));

% What d2 pair is being swapped?
diffmask = squeeze(PolicyPath1(2,:,:,:)) ~= squeeze(PolicyPath2(2,:,:,:));
p1 = squeeze(PolicyPath1(2,:,:,:)); p2 = squeeze(PolicyPath2(2,:,:,:));
fprintf('d2 pairs at diffs: ');
disp(unique([p1(diffmask) p2(diffmask)], 'rows'));


% 
% % Which row of Policy has the diff?
% d = abs(PolicyPath1 - PolicyPath2);
% for r = 1:size(d,1)
%     fprintf('Policy row %d: max abs diff = %g (count = %d)\n', r, ...
%             max(d(r,:),[],'all'), nnz(d(r,:)>0));
% end
% % At first differing state, show both
% mask = any(d>0,1);
% [a_idx,j_idx,t_idx] = ind2sub(size(mask), find(mask,1,'first'));
% fprintf('First diff at (a=%d, j=%d, t=%d, age=%d): nonDC=[%s] DC=[%s]\n', ...
%     a_idx, j_idx, t_idx, j_idx, ...
%     num2str(squeeze(PolicyPath1(:,a_idx,j_idx,t_idx))'), ...
%     num2str(squeeze(PolicyPath2(:,a_idx,j_idx,t_idx))'));
% % And V at the same state (should be ~equal):
% fprintf('V at same state: nonDC=%.15g DC=%.15g\n', ...
%     VPath1(a_idx,j_idx,t_idx), VPath2(a_idx,j_idx,t_idx));


clear VPath1 VPath2 PolicyPath1 PolicyPath2 % PolicyVals1


%% With and without divide-and-conquer (both with fastOLG)
[VPath1,PolicyPath1]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);

% Solve with divide-and-conquer, should give same answer
[VPath2,PolicyPath2]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);

fprintf('Divide-and-conquer (fastOLG), this should be zero: %2.8f \n',max(abs(VPath1(:)-VPath2(:))))
fprintf('Divide-and-conquer (fastOLG), this should be zero: %2.8f \n',max(abs(PolicyPath1(:)-PolicyPath2(:))))

%%
clear VPath1 VPath2 PolicyPath1 PolicyPath2

%% Solve with grid-interpolation. With and without divide-and-conquer (both with slowOLG)
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[VPath3,PolicyPath3]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions3);

% PolicyVals3=PolicyInd2Val_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions3);

% Solve with divide-and-conquer, should give same answer
vfoptions4=vfoptions;
vfoptions4.divideandconquer=1;
vfoptions4.gridinterplayer=1;
vfoptions4.ngridinterp=5;
simoptions4=simoptions;
simoptions4.gridinterplayer=vfoptions4.gridinterplayer;
simoptions4.ngridinterp=vfoptions4.ngridinterp;
[VPath4,PolicyPath4]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions4);

fprintf('Divide-and-conquer (with GI, slowOLG), this should be zero: %2.8f \n',max(abs(VPath3(:)-VPath4(:))))
fprintf('Divide-and-conquer (with GI, slowOLG), this should be zero: %2.8f \n',max(abs(PolicyPath3(:)-PolicyPath4(:))))

clear VPath3 VPath4 PolicyPath3 PolicyPath4 % PolicyVals3


%% Solve with grid-interpolation. With and without divide-and-conquer (both with fastOLG)
[VPath3,PolicyPath3]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions3);

% Solve with divide-and-conquer, should give same answer
[VPath4,PolicyPath4]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final_GI, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);

fprintf('Divide-and-conquer (with GI, fastOLG), this should be zero: %2.8f \n',max(abs(VPath3(:)-VPath4(:))))
fprintf('Divide-and-conquer (with GI, fastOLG), this should be zero: %2.8f \n',max(abs(PolicyPath3(:)-PolicyPath4(:))))

%%
clear VPath3 VPath4 PolicyPath3 PolicyPath4

%% Use a really big a_grid, then the moments should be essentially the same with/without grid interpolation
simoptions2.a_grid=a_grid_big; % ExpAsset codes slice a2_grid out of simoptions.a_grid using n_a1, so must match n_a_big
simoptions4.a_grid=a_grid_big;
[VPath2b,PolicyPath2b]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final_big, Policy_final_big, Params, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
AgentDist_initial_big=StationaryDist_FHorz_Case1(jequaloneDist_big,AgeWeightParamNames,PolicyPath2b(:,:,:,:,1),n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions2);
AgentDistPath2=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_big, jequaloneDist_big, PricePath, ParamPath, PolicyPath2b, AgeWeightParamNames,n_d,n_a_big,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions2);
AggVarsPath2=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath2, PolicyPath2b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, transpathoptionsbaseline, simoptions2);

[VPath4b,PolicyPath4b]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final_big, Policy_final_big_GI, Params, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);
AgentDistPath4=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_big, jequaloneDist_big, PricePath, ParamPath, PolicyPath4b, AgeWeightParamNames,n_d,n_a_big,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions4);
AggVarsPath4=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath4, PolicyPath4b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, N_j, d_grid, a_grid_big,z_grid, transpathoptionsbaseline, simoptions4);
simoptions2.a_grid=a_grid; % restore for the small-grid section below
simoptions4.a_grid=a_grid;

fprintf('With/without grid interp, should get much the same moments (for big a_grid) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %2.8f \n',max(abs(AgentDistPath2(:)-AgentDistPath4(:))))
[AggVarsPath2.earnings.Mean; AggVarsPath4.earnings.Mean]
[AggVarsPath2.assets.Mean; AggVarsPath4.assets.Mean]

%% Do some graphs of the AggVars path to see them
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:T,AggVarsPath2.earnings.Mean, 1:1:T,AggVarsPath4.earnings.Mean)
title('Earnings Mean')
legend('1','2')
subplot(2,1,2); plot(1:1:T,AggVarsPath2.assets.Mean, 1:1:T,AggVarsPath4.assets.Mean)
title('Assets Mean')
legend('1','2')

clear VPath2b VPath4b AggVarsPath4 PolicyPath4b AgentDistPath2 AgentDistPath4

%% If the path is all constant, should just get same answer as when we dont have a TPath
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
AgentDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptions1);
PricePathConstant.r=Params.r*ones(1,T);
ParamPathConstant.sigma=Params.sigma*ones(1,T);
[VPath1,PolicyPath1]=ValueFnOnTransPath_Case1_FHorz(PricePathConstant, ParamPathConstant, T, V1, Policy1, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);
AgentDistPath1=AgentDistOnTransPath_Case1_FHorz(AgentDist1, jequaloneDist, PricePathConstant, ParamPathConstant, PolicyPath1, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions1);
V1_rep=repmat(V1,1,1,1,1,T);
fprintf('Do nothing TPath, this should be zero, V: %2.8f \n',max(abs(VPath1(:)-V1_rep(:))))
Policy1_rep=repmat(Policy1,1,1,1,1,1,T);
fprintf('Do nothing TPath, this should be zero, Policy: %2.8f \n',max(abs(PolicyPath1(:)-Policy1_rep(:))))
AgentDist1_rep=repmat(AgentDist1,1,1,1,1,T);
fprintf('Do nothing TPath, this should be zero, AgentDist: %2.8f \n',max(abs(AgentDistPath1(:)-AgentDist1_rep(:))))

clear V1 Policy1 VPath1 PolicyPath1

%% Run the GE transition path, but with transpathoptions.maxiter=1, so it ends after one iteration -- just a shape-check (the core is tested elsewhere)
[~,PolicyPath2]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyPath2(:,:,:,:,1),n_d,n_a,n_z,N_j,pi_z,Params,simoptions2);
clear PolicyPath2

transpathoptions.maxiter=1;
GeneralEqmEqns.dummy=@(earnings) 0;

transpathoptions.fastOLG=1;
PricePath2=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);

transpathoptions.fastOLG=0;
PricePath2B=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);

fprintf('One iter of TPath, with/without fastOLG, this should be zero: %2.8f \n',max(abs(PricePath2.r-PricePath2B.r)))

% Big grid, uses vfoptions2 with divide-and-conquer
transpathoptions.fastOLG=1;
PricePath3A=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final_big, AgentDist_initial_big, jequaloneDist_big, n_d, n_a_big, n_z, N_j, d_grid,a_grid_big,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions2, vfoptions2);

% vfoptions4 has divide-and-conquer and grid interpolation layer
PricePath3B=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final_big, AgentDist_initial_big, jequaloneDist_big, n_d, n_a_big, n_z, N_j, d_grid,a_grid_big,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions4, vfoptions4);

fprintf('One iter of TPath, with/without GI, this should be close to zero: %2.8f \n',max(abs(PricePath3A.r-PricePath3B.r)))

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
