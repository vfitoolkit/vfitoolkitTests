function output=CoreFHorzPType_NivsZidentity(n_d,n_a,N_j,d_grid,a_grid,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i)
% Test 2: PType with N_i types is equivalent to a non-PType model where z encodes
% the type with an identity transition matrix.
%
% Types differ in the CRRA curvature (sigma):
%   sigma_pt = [2; 3]
% so the per-type parameter dependence enters the ReturnFn directly (not via the
% discount factor, and not via an age-dependent param).
%
% (A) PType solve: N_i=2, no z. Uses ReturnFn_nod_noz_noe_nosemiz with sigma_pt.
% (B) Single solve: n_z=N_i, pi_z=eye(N_i), z_grid=[sigma values]. Uses
%     ReturnFn_zAsSigma which reads z as the curvature.

% PType-specific sigma
Params.sigma_pt=[2; 3];

n_d=0;
d_grid=[];

vfoptions=struct();
simoptions=struct();

%% (A) PType, no z
ReturnFn_A=@(aprime,a,r,w,kappa_j,sigma_pt,agej,Jr,pension) ...
    ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma_pt,agej,Jr,pension);

n_z_A=0;
z_grid_A=[];
pi_z_A=[];

jequaloneDist_A=zeros(n_a,1,'gpuArray');
jequaloneDist_A(1)=1; % no assets

[V_A,Policy_A]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,pi_z_A,ReturnFn_A,Params,DiscountFactorParamNames,vfoptions);
V_A_vfp=ValueFnFromPolicy_FHorz_PType(Policy_A,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,pi_z_A,ReturnFn_A,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_A=StationaryDist_Case1_FHorz_PType(jequaloneDist_A,AgeWeightParamNames,PTypeDistParamNames,Policy_A,n_d,n_a,n_z_A,N_j,N_i,pi_z_A,Params,simoptions);

%% (B) No PType, z encodes the type with identity transitions
n_z_B=N_i;
z_grid_B=Params.sigma_pt; % [2; 3]
pi_z_B=eye(N_i,'gpuArray');

ReturnFn_B=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) ...
    ReturnFn_zAsSigma(aprime,a,z,r,w,kappa_j,agej,Jr,pension);

% Initial dist on (a,z): all mass at a=1, split across types by ptypeweights
jequaloneDist_B=zeros(n_a,n_z_B,'gpuArray');
jequaloneDist_B(1,:)=reshape(Params.(PTypeDistParamNames{:}),1,n_z_B);

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,pi_z_B,ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptions);
V_B_vfp=ValueFnFromPolicy_FHorz(Policy_B,n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,pi_z_B,ReturnFn_B,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a,n_z_B,N_j,pi_z_B,Params,simoptions);

%% Compare per type: V_A.(type ii) ≡ V_B(:,ii,:); StationaryDist_A.(ii)*ptw(ii) ≡ StationaryDist_B(:,ii,:)
names_A=fieldnames(V_A);
ptw=Params.(PTypeDistParamNames{:});
for ii=1:N_i
    nA=names_A{ii};

    V_A_ii=V_A.(nA);                          % [n_a, N_j]
    V_B_ii=reshape(V_B(:,ii,:),[n_a,N_j]);    % strip z-dim

    Pol_A_ii=Policy_A.(nA);
    Pol_B_ii_full=Policy_B(:,:,ii,:);         % strip z-dim (leading 1 dim retained)
    Pol_B_ii=reshape(Pol_B_ii_full,size(Pol_A_ii));

    Dist_A_ii=StationaryDist_A.(nA);
    Dist_B_ii=reshape(StationaryDist_B(:,ii,:),[n_a,N_j]);

    V_A_vfp_ii=V_A_vfp.(nA);
    V_B_vfp_ii=reshape(V_B_vfp(:,ii,:),[n_a,N_j]);

    fprintf('N_i vs z-identity, V    (type %d), this should be zero: %.3e \n',ii,max(abs(V_A_ii(:)-V_B_ii(:))))
    fprintf('N_i vs z-identity, Pol  (type %d), this should be zero: %.3e \n',ii,max(abs(Pol_A_ii(:)-Pol_B_ii(:))))
    fprintf('N_i vs z-identity, Dist (type %d, weighted), this should be zero: %.3e \n',ii,max(abs(ptw(ii)*Dist_A_ii(:)-Dist_B_ii(:))))
    fprintf('N_i vs z-identity, VFP  vs V (PType,   type %d), this should be zero: %.3e \n',ii,max(abs(V_A_vfp_ii(:)-V_A_ii(:))))
    fprintf('N_i vs z-identity, VFP  vs V (noPType, type %d), this should be zero: %.3e \n',ii,max(abs(V_B_vfp_ii(:)-V_B_ii(:))))
end

output=struct();

end
