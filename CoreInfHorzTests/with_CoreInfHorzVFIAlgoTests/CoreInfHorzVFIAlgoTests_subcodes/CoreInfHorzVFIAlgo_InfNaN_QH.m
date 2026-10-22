function output=CoreInfHorzVFIAlgo_InfNaN_QH(n_d,d_grid,n_z,z_grid,pi_z,Params,DiscountFactorParamNames,figure_c)
% Deliberate trigger for the sophisticated quasi-hyperbolic postGI raws (QuasiHyperbolicS_postGI_nod_raw and
% QuasiHyperbolicS_Refine_postGI_raw). Their Howards step iterates Vunderbar at the chosen fine point as
%    (1-wupperH).*EVlowerH+wupperH.*EVupperH
% with EVlowerH/EVupperH the pi-summed (NaN-cleaned) continuations at the two coarse nodes, so a zero weight on a node
% that is -Inf at a reachable z' gives NaN, which lands in Vunderbar and is laundered next pass: the same defect as the
% exponential raws, at the combine after the pi-sum instead of before it.
%
% beta0=1 ONLY (see the engine header): the two raws then run their Howards code exactly as at beta0<1, but the solution
% is exact, so 'vs pure VFI' is an identity. Both devices, nod and d: four figures, figure_c to figure_c+3.

DF=DiscountFactorParamNames;
Params.w=0.5; % so that c<1, u<0, at every state that matters: a continuation priced at 0 is then an inflation of V
Params.beta0=1; % see above
vfoptions=struct();
vfoptions.exoticpreferences='QuasiHyperbolic';
vfoptions.quasi_hyperbolic='Sophisticated';
vfoptions.QHadditionaldiscount='beta0'; % CHAR, as the QH bank notes

%% Device 1: below the natural borrowing limit (fires the zero weight on a dead LOWER node)
n_neg=12;
a_gridNBL=[(-n_neg:1:-1)'; 5*linspace(0,1,101)'.^3];
n_aNBL=length(a_gridNBL);
fprintf('\nbelow NBL, sophisticated quasi-hyperbolic at beta0=1: a_grid has %i integer-spaced nodes from %i to -1 below the usual cubic grid (n_a=%i) \n',n_neg,-n_neg,n_aNBL);

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: nod, z, noe, QH sophisticated beta0=1',0,n_aNBL,n_z,[],a_gridNBL,z_grid,pi_z,ReturnFn,Params,DF,vfoptions,figure_c);

ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: d, z, noe, QH sophisticated beta0=1',n_d,n_aNBL,n_z,d_grid,a_gridNBL,z_grid,pi_z,ReturnFn,Params,DF,vfoptions,figure_c+1);

%% Device 2: a -Inf cap on the state at the top z (fires the zero weight on a dead UPPER node)
Params.acap=0;
Params.zcap=(z_grid(end-1)+z_grid(end))/2;
a_gridcap=5*linspace(0,1,101)'.^3;
n_acap=length(a_gridcap);
pi_zcap=pi_z;
pi_zcap(1:end-2,end)=0; % the three lowest z never reach z_top
pi_zcap=pi_zcap./sum(pi_zcap,2);
fprintf('\ncap at top z, sophisticated quasi-hyperbolic at beta0=1: V(a,z_top)=-Inf for a>%g; pi_z column %i zeroed in rows 1:%i and rows renormalised \n',Params.acap,n_z,n_z-2);

ReturnFn=@(aprime,a,z,r,w,sigma,acap,zcap) ReturnFn_nod_z_noe_InfNaNcap(aprime,a,z,r,w,sigma,acap,zcap);
output=CoreInfHorzVFIAlgo_InfNaNcompare('cap at top z: nod, z, noe, QH sophisticated beta0=1',0,n_acap,n_z,[],a_gridcap,z_grid,pi_zcap,ReturnFn,Params,DF,vfoptions,figure_c+2);

ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi,acap,zcap) ReturnFn_d_z_noe_InfNaNcap(d,aprime,a,z,r,w,sigma,eta,varphi,acap,zcap);
output=CoreInfHorzVFIAlgo_InfNaNcompare('cap at top z: d, z, noe, QH sophisticated beta0=1',n_d,n_acap,n_z,d_grid,a_gridcap,z_grid,pi_zcap,ReturnFn,Params,DF,vfoptions,figure_c+3);

end
