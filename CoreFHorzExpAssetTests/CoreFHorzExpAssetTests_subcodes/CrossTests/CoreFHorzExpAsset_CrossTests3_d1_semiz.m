function output=CoreFHorzExpAsset_CrossTests3_d1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% Just solve without z and without e
ReturnFn_twoendo=@(d1,d3,a1prime,a2prime,a1,a2,semiz,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost)...
    ReturnFn_TwoEndo_d1_noz_noe_semiz(d1,d3,a1prime,a2prime,a1,a2,semiz,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);

ReturnFn_none=@(d1,d2,d3,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost)...
    ReturnFn_d1_noz_noe_semiz(d1,d2,d3,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);

% Note: d2 and a2prime should end up just the same thing. The two endo
% problem is kind of stupid, but that it okay as it is just to check the
% compute.

aprimeFn=@(d2,a2) d2; % a2prime is just d2
d_grid=[d_grid(1:n_d(1)); a_grid(n_a(1)+1:end); d_grid(n_d(1)+n_d(2)+1:end)]; % keep d1&d3, set d2 to a2prime (a2prime is just a2)
n_d(2)=n_a(2); % keep d1&d3, set d2 to a2prime (a2prime is just a2)

%% Solving with just a single points for z with value 1 and prob 1 gives us same as no shocks
jequaloneDist_none=zeros([n_a,simoptionsbaseline.n_semiz],'gpuArray');
jequaloneDist_none(1,1,ceil(simoptionsbaseline.n_semiz/2))=1; % no assets

n_d_alt=[n_d(1),n_d(3)]; % keep d1&d3, drop d2
d_grid_alt=[d_grid(1:n_d(1)); d_grid(n_d(1)+n_d(2)+1:end)];  % keep d1&d3, drop d2

vfoptionsA.divideandconquer=0;
% Setup semiz
vfoptionsA.n_semiz=vfoptionsbaseline.n_semiz;
vfoptionsA.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptionsA.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptionsA.n_semiz=simoptionsbaseline.n_semiz;
simoptionsA.semiz_grid=simoptionsbaseline.semiz_grid;
simoptionsA.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptionsA.d_grid=d_grid_alt;
[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d_alt,n_a,0,N_j,d_grid_alt,a_grid,[],[],ReturnFn_twoendo,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d_alt,n_a,0,N_j,[],Params,simoptionsA);

Policy0alt=[Policy0(1,:,:,:); Policy0(4,:,:,:); Policy0(2,:,:,:); Policy0(3,:,:,:)]; % d1, a2prime (which is d2), d3, a1prime, 

% Experience asset
vfoptionsB.experienceasset=1;
simoptionsB.experienceasset=1;
vfoptionsB.aprimeFn=aprimeFn;
simoptionsB.aprimeFn=aprimeFn;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;
% Setup semiz
vfoptionsB.n_semiz=vfoptionsbaseline.n_semiz;
vfoptionsB.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptionsB.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptionsB.n_semiz=simoptionsbaseline.n_semiz;
simoptionsB.semiz_grid=simoptionsbaseline.semiz_grid;
simoptionsB.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptionsB.d_grid=d_grid;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy1,n_d,n_a,0,N_j,1,Params,simoptionsB);

fprintf('Cross test 3: expasset is just a standard endo state, this should be zero: %.3e \n',max(abs(V0(:)-V1(:))))
fprintf('Cross test 3: expasset is just a standard endo state, this should be zero: %.3e \n',max(abs(Policy0alt(:)-Policy1(:))))
fprintf('Cross test 3: expasset is just a standard endo state, this should be zero: %.3e \n',max(abs(StationaryDist0(:)-StationaryDist1(:))))


%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end