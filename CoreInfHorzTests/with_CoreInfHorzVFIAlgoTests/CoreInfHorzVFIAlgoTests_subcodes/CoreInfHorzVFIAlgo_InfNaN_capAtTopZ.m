function output=CoreInfHorzVFIAlgo_InfNaN_capAtTopZ(n_d,d_grid,n_z,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,figure_c)
% Deliberate trigger, device 2: a dead region ABOVE finite states, at one z only.
%
% The ReturnFn returns -Inf for every a above acap=a_grid(1)=0 at the top z state, whatever a' is
% (ReturnFn_*_InfNaNcap: a cap on the STATE, not on a'), so V(a,z_top)=-Inf for a>0 while V(a,z) at
% the other z is finite. A state that can reach z_top can only choose a'=0 (any a'>0 risks a -Inf),
% so its chosen point is the bottom node of the grid, which is the bottom node of its +-maxaprimediff
% window: the one place the postGI raws code with all the weight on the LOWER node (lowerprob==1).
% The Howards step then puts a ZERO weight on the node above, which is -Inf at z'=z_top:
% 0*(-Inf)=NaN, and the (a,z,z_top) continuation is priced at 0. This device fires the 'dead upper
% node' branch of the defect (the EVupper(lowerprob==1)=0 line of the fix).
%
% So that the model is not everyone-at-zero: the three lowest z states never move to z_top (that
% column of pi_z is zeroed and the rows renormalised), so they save freely; the top two z states are
% the ones forced to a'=0. (The proposal wrote the cap as a'<=abar; that leaves V finite above abar,
% since the agent just consumes down, so it is on a that the cap has to go.)
%
% The without-z models have no z to cap at, and are not run here (device 1 covers their raws).
% Six figures, figure_c to figure_c+5.

DF=DiscountFactorParamNames;
Params.w=0.5; % so that c<1, u<0, at the forced states: a continuation priced at 0 is then an inflation of V
Params.acap=0;
Params.zcap=(z_grid(end-1)+z_grid(end))/2; % between the top two z: the cap bites at z_top only

a_grid=5*linspace(0,1,101)'.^3; % the bank's usual grid; a_grid(1)=0=acap is the only alive node at z_top
n_a=length(a_grid);

pi_z(1:end-2,end)=0; % the three lowest z never reach z_top
pi_z=pi_z./sum(pi_z,2);
fprintf('\ncap at top z: V(a,z_top)=-Inf for a>%g; pi_z column %i zeroed in rows 1:%i and rows renormalised \n',Params.acap,n_z,n_z-2);

vfoptions=struct();
vfoptions_e=struct();
vfoptions_e.n_e=vfoptionsbaseline.n_e;
vfoptions_e.e_grid=vfoptionsbaseline.e_grid;
vfoptions_e.pi_e=vfoptionsbaseline.pi_e;

%% nod, z, noe (postGI_nod_raw, postGI_sparse_nod_raw, HowardMix_nod_raw)
ReturnFn=@(aprime,a,z,r,w,sigma,acap,zcap) ReturnFn_nod_z_noe_InfNaNcap(aprime,a,z,r,w,sigma,acap,zcap);
output=CoreInfHorzVFIAlgo_InfNaNcompare('cap at top z: nod, z, noe',0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions,figure_c);

%% d, z, noe (Refine_postGI_raw, postGI_sparse_raw, Refine_postGI_HowardMix_raw)
ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi,acap,zcap) ReturnFn_d_z_noe_InfNaNcap(d,aprime,a,z,r,w,sigma,eta,varphi,acap,zcap);
output=CoreInfHorzVFIAlgo_InfNaNcompare('cap at top z: d, z, noe',n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions,figure_c+1);

%% nod, z, e (POSITIVE CONTROL: the _e raws already zero each weighted term)
ReturnFn=@(aprime,a,z,e,r,w,sigma,acap,zcap) ReturnFn_nod_z_e_InfNaNcap(aprime,a,z,e,r,w,sigma,acap,zcap);
output=CoreInfHorzVFIAlgo_InfNaNcompare('cap at top z: nod, z, e (positive control, the _e raws have the per-term form)',0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions_e,figure_c+2);

%% d, z, e (POSITIVE CONTROL)
ReturnFn=@(d,aprime,a,z,e,r,w,sigma,eta,varphi,acap,zcap) ReturnFn_d_z_e_InfNaNcap(d,aprime,a,z,e,r,w,sigma,eta,varphi,acap,zcap);
output=CoreInfHorzVFIAlgo_InfNaNcompare('cap at top z: d, z, e (positive control, the _e raws have the per-term form)',n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions_e,figure_c+3);

%% Two endogenous states: the cap is on a1 (the grid interpolation layer is on a1); a2 is free, so the
% forced states still have a real savings choice in a2
n_a2A=[31,11];
a1_grid2A=5*linspace(0,1,31)'.^3;
a2_grid2A=3*linspace(0,1,11)'.^3;
a_grid2A=[a1_grid2A; a2_grid2A]; % stacked, as the toolkit expects for multiple endogenous states
Params.r2=0.03; % return on the second asset (the first pays Params.r)

%% 2A nod, z (postGI2A_nod_raw, postGI2A_sparse_nod_raw)
ReturnFn=@(a1prime,a2prime,a1,a2,z,r,r2,w,sigma,acap,zcap) ReturnFn_nod_z_noe_InfNaNcap_with2A(a1prime,a2prime,a1,a2,z,r,r2,w,sigma,acap,zcap);
output=CoreInfHorzVFIAlgo_InfNaNcompare('cap at top z: 2A nod, z, noe',0,n_a2A,n_z,[],a_grid2A,z_grid,pi_z,ReturnFn,Params,DF,vfoptions,figure_c+4);

%% 2A d, z (Refine_postGI2A_raw, Refine_postGI2A_sparse_raw)
ReturnFn=@(d,a1prime,a2prime,a1,a2,z,r,r2,w,sigma,eta,varphi,acap,zcap) ReturnFn_d_z_noe_InfNaNcap_with2A(d,a1prime,a2prime,a1,a2,z,r,r2,w,sigma,eta,varphi,acap,zcap);
output=CoreInfHorzVFIAlgo_InfNaNcompare('cap at top z: 2A d, z, noe',n_d,n_a2A,n_z,d_grid,a_grid2A,z_grid,pi_z,ReturnFn,Params,DF,vfoptions,figure_c+5);

end
