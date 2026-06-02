function output=CoreFHorzPType_ShockTests_8types(n_d,n_a,n_z,n_d_semiz,d_grid_semiz,n_d2_semiz,d2_grid_semiz,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,vfoptionsbaseline,simoptionsbaseline)
% Shock-combination cross test for PType. One PType solve with 8 types covering
% every (z, e, semiz) combination should give the same per-type V, Policy and
% StationaryDist as 8 separate non-PType solves.
%
% Per-type n_d, d_grid, n_z, z_grid, pi_z, ReturnFn, jequaloneDist and the
% e/semiz pieces of vfoptions/simoptions are all passed as Names_i-keyed
% structures (required when these differ across types).
%
% All semiz types use the same nod1 + binary d2 layout used in CoreFHorzTests.

Names_i={'none','z','e','semiz','ze','zsemiz','esemiz','zesemiz'};
N_i=length(Names_i);

% Arbitrary unequal type weights
Params.ptypeweights=(1:N_i)'/sum(1:N_i);

% Pull the e and semiz pieces out of vfoptionsbaseline for convenience
n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;
n_semiz=vfoptionsbaseline.n_semiz;
semiz_grid=vfoptionsbaseline.semiz_grid;
SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

%% Per-type structs for the single PType call
n_d_PT=struct();           d_grid_PT=struct();
n_z_PT=struct();           z_grid_PT=struct();           pi_z_PT=struct();
ReturnFn_PT=struct();      jequaloneDist_PT=struct();

vfopt_PT=struct();         simopt_PT=struct();
vfopt_PT.n_e=struct();     vfopt_PT.e_grid=struct();     vfopt_PT.pi_e=struct();
vfopt_PT.n_semiz=struct(); vfopt_PT.semiz_grid=struct(); vfopt_PT.SemiExoStateFn=struct();
simopt_PT.n_e=struct();    simopt_PT.e_grid=struct();    simopt_PT.pi_e=struct();
simopt_PT.n_semiz=struct();simopt_PT.semiz_grid=struct();simopt_PT.SemiExoStateFn=struct();

% type 1: none
nm='none';
n_d_PT.(nm)=0; d_grid_PT.(nm)=[];
n_z_PT.(nm)=0; z_grid_PT.(nm)=[]; pi_z_PT.(nm)=[];
ReturnFn_PT.(nm)=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);
jequaloneDist_PT.(nm)=zeros(n_a,1,'gpuArray');
jequaloneDist_PT.(nm)(1)=1;

% type 2: z only
nm='z';
n_d_PT.(nm)=0; d_grid_PT.(nm)=[];
n_z_PT.(nm)=n_z; z_grid_PT.(nm)=z_grid; pi_z_PT.(nm)=pi_z;
ReturnFn_PT.(nm)=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
jequaloneDist_PT.(nm)=zeros(n_a,n_z,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_z/2))=1;

% type 3: e only
nm='e';
n_d_PT.(nm)=0; d_grid_PT.(nm)=[];
n_z_PT.(nm)=0; z_grid_PT.(nm)=[]; pi_z_PT.(nm)=[];
ReturnFn_PT.(nm)=@(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension);
vfopt_PT.n_e.(nm)=n_e; vfopt_PT.e_grid.(nm)=e_grid; vfopt_PT.pi_e.(nm)=pi_e;
simopt_PT.n_e.(nm)=n_e;simopt_PT.e_grid.(nm)=e_grid;simopt_PT.pi_e.(nm)=pi_e;
jequaloneDist_PT.(nm)=zeros(n_a,n_e,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_e/2))=1;

% type 4: semiz only
nm='semiz';
n_d_PT.(nm)=n_d2_semiz; d_grid_PT.(nm)=d2_grid_semiz;
n_z_PT.(nm)=0; z_grid_PT.(nm)=[]; pi_z_PT.(nm)=[];
ReturnFn_PT.(nm)=@(d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ...
    ReturnFn_nod1_noz_noe_semiz(d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
vfopt_PT.n_semiz.(nm)=n_semiz; vfopt_PT.semiz_grid.(nm)=semiz_grid; vfopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
simopt_PT.n_semiz.(nm)=n_semiz;simopt_PT.semiz_grid.(nm)=semiz_grid;simopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
jequaloneDist_PT.(nm)=zeros(n_a,n_semiz,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_semiz/2))=1;

% type 5: z + e
nm='ze';
n_d_PT.(nm)=0; d_grid_PT.(nm)=[];
n_z_PT.(nm)=n_z; z_grid_PT.(nm)=z_grid; pi_z_PT.(nm)=pi_z;
ReturnFn_PT.(nm)=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
vfopt_PT.n_e.(nm)=n_e; vfopt_PT.e_grid.(nm)=e_grid; vfopt_PT.pi_e.(nm)=pi_e;
simopt_PT.n_e.(nm)=n_e;simopt_PT.e_grid.(nm)=e_grid;simopt_PT.pi_e.(nm)=pi_e;
jequaloneDist_PT.(nm)=zeros(n_a,n_z,n_e,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_z/2),ceil(n_e/2))=1;

% type 6: z + semiz
nm='zsemiz';
n_d_PT.(nm)=n_d2_semiz; d_grid_PT.(nm)=d2_grid_semiz;
n_z_PT.(nm)=n_z; z_grid_PT.(nm)=z_grid; pi_z_PT.(nm)=pi_z;
ReturnFn_PT.(nm)=@(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ...
    ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
vfopt_PT.n_semiz.(nm)=n_semiz; vfopt_PT.semiz_grid.(nm)=semiz_grid; vfopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
simopt_PT.n_semiz.(nm)=n_semiz;simopt_PT.semiz_grid.(nm)=semiz_grid;simopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
jequaloneDist_PT.(nm)=zeros(n_a,n_semiz,n_z,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_semiz/2),ceil(n_z/2))=1;

% type 7: e + semiz
nm='esemiz';
n_d_PT.(nm)=n_d2_semiz; d_grid_PT.(nm)=d2_grid_semiz;
n_z_PT.(nm)=0; z_grid_PT.(nm)=[]; pi_z_PT.(nm)=[];
ReturnFn_PT.(nm)=@(d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ...
    ReturnFn_nod1_noz_e_semiz(d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
vfopt_PT.n_e.(nm)=n_e; vfopt_PT.e_grid.(nm)=e_grid; vfopt_PT.pi_e.(nm)=pi_e;
simopt_PT.n_e.(nm)=n_e;simopt_PT.e_grid.(nm)=e_grid;simopt_PT.pi_e.(nm)=pi_e;
vfopt_PT.n_semiz.(nm)=n_semiz; vfopt_PT.semiz_grid.(nm)=semiz_grid; vfopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
simopt_PT.n_semiz.(nm)=n_semiz;simopt_PT.semiz_grid.(nm)=semiz_grid;simopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
jequaloneDist_PT.(nm)=zeros(n_a,n_semiz,n_e,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_semiz/2),ceil(n_e/2))=1;

% type 8: z + e + semiz
nm='zesemiz';
n_d_PT.(nm)=n_d2_semiz; d_grid_PT.(nm)=d2_grid_semiz;
n_z_PT.(nm)=n_z; z_grid_PT.(nm)=z_grid; pi_z_PT.(nm)=pi_z;
ReturnFn_PT.(nm)=@(d2,aprime,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ...
    ReturnFn_nod1_z_e_semiz(d2,aprime,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
vfopt_PT.n_e.(nm)=n_e; vfopt_PT.e_grid.(nm)=e_grid; vfopt_PT.pi_e.(nm)=pi_e;
simopt_PT.n_e.(nm)=n_e;simopt_PT.e_grid.(nm)=e_grid;simopt_PT.pi_e.(nm)=pi_e;
vfopt_PT.n_semiz.(nm)=n_semiz; vfopt_PT.semiz_grid.(nm)=semiz_grid; vfopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
simopt_PT.n_semiz.(nm)=n_semiz;simopt_PT.semiz_grid.(nm)=semiz_grid;simopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
jequaloneDist_PT.(nm)=zeros(n_a,n_semiz,n_z,n_e,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_semiz/2),ceil(n_z/2),ceil(n_e/2))=1;

%% PType solve (one call, eight types)
[V_PT,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d_PT,n_a,n_z_PT,N_j,Names_i,d_grid_PT,a_grid,z_grid_PT,pi_z_PT,ReturnFn_PT,Params,DiscountFactorParamNames,vfopt_PT);
simopt_PT.d_grid=d_grid_PT;
StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist_PT,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d_PT,n_a,n_z_PT,N_j,Names_i,pi_z_PT,Params,simopt_PT);

%% Eight independent single-type solves
V_solo=struct(); Policy_solo=struct(); Dist_solo=struct();
for ii=1:N_i
    nm=Names_i{ii};

    % Per-type vfoptions/simoptions for the non-PType call
    vfopt_solo=struct();
    simopt_solo=struct();
    if isfield(vfopt_PT.n_e,nm)
        vfopt_solo.n_e=vfopt_PT.n_e.(nm);
        vfopt_solo.e_grid=vfopt_PT.e_grid.(nm);
        vfopt_solo.pi_e=vfopt_PT.pi_e.(nm);
        simopt_solo.n_e=vfopt_solo.n_e;
        simopt_solo.e_grid=vfopt_solo.e_grid;
        simopt_solo.pi_e=vfopt_solo.pi_e;
    end
    if isfield(vfopt_PT.n_semiz,nm)
        vfopt_solo.n_semiz=vfopt_PT.n_semiz.(nm);
        vfopt_solo.semiz_grid=vfopt_PT.semiz_grid.(nm);
        vfopt_solo.SemiExoStateFn=vfopt_PT.SemiExoStateFn.(nm);
        simopt_solo.n_semiz=vfopt_solo.n_semiz;
        simopt_solo.semiz_grid=vfopt_solo.semiz_grid;
        simopt_solo.SemiExoStateFn=vfopt_solo.SemiExoStateFn;
        simopt_solo.d_grid=d_grid_PT.(nm);
    end

    [V_solo.(nm),Policy_solo.(nm)]=ValueFnIter_Case1_FHorz(n_d_PT.(nm),n_a,n_z_PT.(nm),N_j,d_grid_PT.(nm),a_grid,z_grid_PT.(nm),pi_z_PT.(nm),ReturnFn_PT.(nm),Params,DiscountFactorParamNames,[],vfopt_solo);
    Dist_solo.(nm)=StationaryDist_FHorz_Case1(jequaloneDist_PT.(nm),AgeWeightParamNames,Policy_solo.(nm),n_d_PT.(nm),n_a,n_z_PT.(nm),N_j,pi_z_PT.(nm),Params,simopt_solo);
end

%% Compare per type
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('ShockTests 8 types, V    (type %s), this should be zero: %2.8f \n',nm,max(abs(V_PT.(nm)(:)-V_solo.(nm)(:))))
    fprintf('ShockTests 8 types, Pol  (type %s), this should be zero: %2.8f \n',nm,max(abs(Policy_PT.(nm)(:)-Policy_solo.(nm)(:))))
    fprintf('ShockTests 8 types, Dist (type %s), this should be zero: %2.8f \n',nm,max(abs(StationaryDist_PT.(nm)(:)-Dist_solo.(nm)(:))))
end

output=struct();

end
