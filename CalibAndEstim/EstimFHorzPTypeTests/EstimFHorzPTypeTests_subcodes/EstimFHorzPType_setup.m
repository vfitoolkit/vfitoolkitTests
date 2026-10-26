% The model shared by EstimFHorzPTypeTests and its sub-bank withMaxLikelihood/ (moved out of EstimFHorzPTypeTests.m 2026-10-16, unchanged).
% A script, not a function: it creates the variables in the caller's workspace.
%% The model: the CoreFHorzPTypeTests part-22 model on small grids
n_d=0;
d_grid=[];
n_a=101;
a_grid=5*linspace(0,1,n_a)'.^3;
n_z=5;
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.03,n_z);
z_grid=exp(z_grid);
N_j=20;
Names_i={'low','high'};

Params.beta=0.95; % discount factor
DiscountFactorParamNames={'beta'};
Params.mewj=0.9.^(0:N_j-1); Params.mewj=Params.mewj/sum(Params.mewj);
AgeWeightParamNames={'mewj'};
Params.ptypeweights=[0.7;0.3];
PTypeDistParamNames={'ptypeweights'};
Params.sigma=2; % CES utility param for consumption
Params.w=1;
Params.r=0.05;
Params.Jr=16;
Params.pension=0.5;
Params.agej=1:1:N_j;
Params.kappa_j=[0.5:0.1:1,ones(1,9),zeros(1,5)]; % deterministic earnings profile, zero in retirement
Params.kappa_j_pt=[Params.kappa_j; 1.5*Params.kappa_j]; % per type: N_i is the first dimension

ReturnFn=@(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension);

jequaloneDist=zeros(n_a,n_z,'gpuArray'); % everyone starts with no assets at the middle shock (both types)
jequaloneDist(1,ceil(n_z/2))=1;
