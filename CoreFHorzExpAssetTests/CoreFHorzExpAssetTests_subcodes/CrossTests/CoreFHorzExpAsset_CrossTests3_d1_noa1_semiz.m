function output=CoreFHorzExpAsset_CrossTests3_d1_noa1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test 3 noa1+semiz d1: ExpAsset noa1+d1+semiz with degenerate aprimeFn(d2,a2)=d2 should match
% standard 1-endo Case1 with d1+semiz (where aprime plays the role of d2).
% n_d=[n_d1, n_d2, n_d3]; d_grid=[d1_grid; d2_grid; d3_grid].

n_d1_only=n_d(1);
d1_grid_only=d_grid(1:n_d(1));
n_d3_only=n_d(3);
d3_grid_only=d_grid(n_d(1)+n_d(2)+1:end);

semizopts=struct();
semizopts.n_semiz=vfoptionsbaseline.n_semiz;
semizopts.semiz_grid=vfoptionsbaseline.semiz_grid;
semizopts.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

ReturnFn_oneendo=@(d1,d3,aprime,a,semiz,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_OneEndo_d1_noz_noe_noa1_semiz(d1,d3,aprime,a,semiz,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_expasset=@(d1,d2,d3,a,semiz,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_d1_noz_noe_noa1_semiz(d1,d2,d3,a,semiz,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);

aprimeFn=@(d2,a2) d2;

% Side A: standard Case1 with d=[d1,d3]; aprime is the additional choice over a_grid
n_d_oneendo=[n_d1_only, n_d3_only];
d_grid_oneendo=[d1_grid_only; d3_grid_only];

% Side B: ExpAsset noa1; d = [d1, d2(=aprime grid), d3]
n_d_expasset=[n_d1_only, n_a, n_d3_only];
d_grid_expasset=[d1_grid_only; a_grid; d3_grid_only];

jequaloneDist=zeros([n_a,semizopts.n_semiz],'gpuArray');
jequaloneDist(1,ceil(semizopts.n_semiz/2))=1;

% Side A
vfoptionsA=semizopts;
simoptionsA=semizopts;
simoptionsA.d_grid=d_grid_oneendo;
simoptionsA.a_grid=a_grid;
[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d_oneendo,n_a,0,N_j,d_grid_oneendo,a_grid,[],[],ReturnFn_oneendo,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy0,n_d_oneendo,n_a,0,N_j,[],Params,simoptionsA);

% Side B
vfoptionsB=semizopts;
vfoptionsB.experienceasset=1;
vfoptionsB.aprimeFn=aprimeFn;
simoptionsB=semizopts;
simoptionsB.experienceasset=1;
simoptionsB.aprimeFn=aprimeFn;
simoptionsB.d_grid=d_grid_expasset;
simoptionsB.a_grid=a_grid;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d_expasset,n_a,0,N_j,d_grid_expasset,a_grid,[],[],ReturnFn_expasset,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1,n_d_expasset,n_a,0,N_j,1,Params,simoptionsB);

% Side A stores Policy rows as (d1, d3, aprime); Side B (expasset) stores them as (d1, d2, d3),
% where d2 is the a_grid choice that plays aprime's role via aprimeFn=@(d2,a2)d2. The two
% formulations thus put the same decisions in different rows, so reorder Side B to
% (d1, d3, d2=aprime) before comparing -- otherwise the raw element-wise diff is n_a-1 (=12
% here), from lining d3 (1..2) up against aprime (1..n_a). (V and Dist already confirm the two
% models coincide.)
Policy1_aligned=Policy1([1 3 2],:,:,:); % (d1,d2,d3) -> (d1,d3,d2=aprime), matching Policy0's (d1,d3,aprime)
fprintf('Cross test 3 (noa1+d1+semiz): expasset noa1+d1+semiz is just a standard 1-endo+d1+semiz state, this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',max(abs(V0(:)-V1(:))),max(abs(Policy0(:)-Policy1_aligned(:))),max(abs(StationaryDist0(:)-StationaryDist1(:))))

output=struct();

end
