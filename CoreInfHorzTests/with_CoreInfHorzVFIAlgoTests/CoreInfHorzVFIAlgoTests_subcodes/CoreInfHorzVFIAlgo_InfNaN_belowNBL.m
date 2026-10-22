function output=CoreInfHorzVFIAlgo_InfNaN_belowNBL(n_d,d_grid,n_z,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,figure_c)
% Deliberate trigger, device 1: the asset grid extends BELOW the natural borrowing limit.
%
% With w=0.5 and r=0.05 the natural borrowing limit at the lowest z (0.85) is -w*z/r=-8.5, so on an
% integer-spaced block -12,...,-1 below the usual cubic grid on [0,5] the nodes -12 to -9 are dead
% ends at every z (no a' keeps c>0 and V(a',z') finite at every reachable z'), while -8 is alive:
% at a=-8, a'=-8 and z=0.85, c=0.025. (With a labour choice d<1 the top of d_grid is 0.875, so the
% limit moves to -7.4 and -8 is dead too; with e the lowest e is 0.84, limit -7.1; without z the
% limit is -10. The no-GI solve defines the dead set, the engine reports it.)
%
% The constrained states, the poor ones who would borrow more, choose the first alive node EXACTLY,
% because every fine point between the dead node and the alive node interpolates a -Inf. The postGI
% raws code a chosen coarse node as (lower=node-1, all weight on the upper node), so the Howards
% step puts a ZERO weight on the dead node: 0*(-Inf)=NaN, the whole (a,z,z') entry is NaN, and the
% isnan line prices that continuation at 0. This device fires the 'dead lower node' branch of the
% defect (the EVlower(lowerprob==0)=0 line of the fix), at every z row where the constraint binds.
%
% Runs every model variant the eight base postGI raws need: nod/d, with/without z, with e (the _e
% raws already have the per-term form, so those two are a POSITIVE CONTROL: they must be green even
% before the fix), and two endogenous states (the GI2A raws). Eight figures, figure_c to figure_c+7.
%
% Uses the bank's own ReturnFns: the device is the grid, not the return function.

DF=DiscountFactorParamNames;
Params.w=0.5; % so that c<1, u<0, at every state that matters: a continuation priced at 0 is then an inflation of V

n_neg=12;
a_grid=[(-n_neg:1:-1)'; 5*linspace(0,1,101)'.^3];
n_a=length(a_grid);
fprintf('\nbelow NBL: a_grid has %i integer-spaced nodes from %i to -1 below the usual cubic grid on [0,5] (n_a=%i); natural borrowing limit at the lowest z, w=0.5, r=0.05: %.2f \n',n_neg,-n_neg,n_a,-Params.w*min(z_grid)/Params.r);

vfoptions=struct();
vfoptions_e=struct();
vfoptions_e.n_e=vfoptionsbaseline.n_e;
vfoptions_e.e_grid=vfoptionsbaseline.e_grid;
vfoptions_e.pi_e=vfoptionsbaseline.pi_e;

%% nod, z, noe (postGI_nod_raw, postGI_sparse_nod_raw, HowardMix_nod_raw)
ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: nod, z, noe',0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions,figure_c);

%% d, z, noe (Refine_postGI_raw, postGI_sparse_raw, Refine_postGI_HowardMix_raw)
ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: d, z, noe',n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions,figure_c+1);

%% nod, z, e (POSITIVE CONTROL: postGI_nod_e_raw and postGI_sparse_nod_e_raw already zero each weighted term)
ReturnFn=@(aprime,a,z,e,r,w,sigma) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,sigma);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: nod, z, e (positive control, the _e raws have the per-term form)',0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions_e,figure_c+2);

%% d, z, e (POSITIVE CONTROL: Refine_postGI_e_raw and postGI_sparse_e_raw)
ReturnFn=@(d,aprime,a,z,e,r,w,sigma,eta,varphi) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,sigma,eta,varphi);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: d, z, e (positive control, the _e raws have the per-term form)',n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions_e,figure_c+3);

%% nod, noz, noe (postGI_nod_noz_raw)
ReturnFn=@(aprime,a,r,w,sigma) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,sigma);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: nod, noz, noe',0,n_a,0,[],a_grid,[],[],ReturnFn,Params,DF,vfoptions,figure_c+4);

%% d, noz, noe (Refine_postGI_noz_raw)
ReturnFn=@(d,aprime,a,r,w,sigma,eta,varphi) ReturnFn_d_noz_noe_nosemiz(d,aprime,a,r,w,sigma,eta,varphi);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: d, noz, noe',n_d,n_a,0,d_grid,a_grid,[],[],ReturnFn,Params,DF,vfoptions,figure_c+5);

%% Two endogenous states: the grid interpolation layer is on a1, so the dead region sits on the a1 grid
% (as in Part 3b of algocompare: small grids, this is a correctness check not a timing one)
n_a2A=[n_neg+31,11];
a1_grid2A=[(-n_neg:1:-1)'; 5*linspace(0,1,31)'.^3];
a2_grid2A=3*linspace(0,1,11)'.^3;
a_grid2A=[a1_grid2A; a2_grid2A]; % stacked, as the toolkit expects for multiple endogenous states
Params.r2=0.03; % return on the second asset (the first pays Params.r)

%% 2A nod, z (postGI2A_nod_raw, postGI2A_sparse_nod_raw)
ReturnFn=@(a1prime,a2prime,a1,a2,z,r,r2,w,sigma) ReturnFn_nod_z_noe_nosemiz_with2A(a1prime,a2prime,a1,a2,z,r,r2,w,sigma);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: 2A nod, z, noe',0,n_a2A,n_z,[],a_grid2A,z_grid,pi_z,ReturnFn,Params,DF,vfoptions,figure_c+6);

%% 2A d, z (Refine_postGI2A_raw, Refine_postGI2A_sparse_raw)
ReturnFn=@(d,a1prime,a2prime,a1,a2,z,r,r2,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz_with2A(d,a1prime,a2prime,a1,a2,z,r,r2,w,sigma,eta,varphi);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: 2A d, z, noe',n_d,n_a2A,n_z,d_grid,a_grid2A,z_grid,pi_z,ReturnFn,Params,DF,vfoptions,figure_c+7);

end
