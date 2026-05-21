function output=CoreFHorzTwoEndo_CrossTests_d1_semiz(n_d_semiz,n_d2_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross tests for the two-endo SemiExo (d1) code paths.
% Test 4: with-d1 model where d1 takes 3 values but the ReturnFn ignores it.
% V from the with-d1 code path must equal V from the nod1 code path.

%% Test 4: d1 ignored in ReturnFn -> with-d1 path equals nod1 path
n_d1_t6=3;
d1_grid_t6=[0.1; 0.5; 0.9]; % values don't matter; ReturnFn ignores d1
n_d_t6=[n_d1_t6, n_d2_semiz(end)]; % [n_d1, n_d2]
d_grid_t6=[d1_grid_t6; d2_grid_semiz]; % stacked grid [d1_grid; d2_grid]

vfopts=struct();
vfopts.n_semiz=vfoptionsbaseline.n_semiz;
vfopts.semiz_grid=vfoptionsbaseline.semiz_grid;
vfopts.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

% with-d1 ReturnFn (d1 in signature but ignored)
ReturnFn_d1=@(d1,d2,a1prime,a2prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2) ReturnFn_d1_noz_noe_semiz_test6(d1,d2,a1prime,a2prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2);

[V_d1,Pol_d1]=ValueFnIter_Case1_FHorz(n_d_t6,n_a,0,N_j,d_grid_t6,a_grid,[],[],ReturnFn_d1,Params,DiscountFactorParamNames,[],vfopts);

% nod1 reference: only d2, standard nod1 semiz ReturnFn
ReturnFn_nod1=@(d2,a1prime,a2prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2) ReturnFn_nod1_noz_noe_semiz(d2,a1prime,a2prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2);
[V_nod1,Pol_nod1]=ValueFnIter_Case1_FHorz(n_d2_semiz,n_a,0,N_j,d2_grid_semiz,a_grid,[],[],ReturnFn_nod1,Params,DiscountFactorParamNames,[],vfopts);

fprintf('Test 4 (d1 ignored vs nod1), V should be zero: %2.8f \n',max(abs(V_d1(:)-V_nod1(:))))

%%
output=struct(); % not currently used

end
