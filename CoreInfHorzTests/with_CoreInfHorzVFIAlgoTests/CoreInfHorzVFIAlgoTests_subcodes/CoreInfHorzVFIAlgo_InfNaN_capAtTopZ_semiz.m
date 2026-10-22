function output=CoreInfHorzVFIAlgo_InfNaN_capAtTopZ_semiz(n_d2_semiz,n_d_semiz,d2_grid_semiz,d_grid_semiz,n_z,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,figure_c)
% Deliberate trigger, device 2 (a -Inf cap on the state at the top z, so the forced choice a'=0 is the window's bottom
% node and the zero weight lands on a dead UPPER node), for the semiz models with z. See CoreInfHorzVFIAlgo_InfNaN_capAtTopZ
% for the device. Four variants: nod1/d1, e/noe. Four figures, figure_c to figure_c+3.

DF=DiscountFactorParamNames;
Params.w=0.5; % so that c<1, u<0, at the forced states: a continuation priced at 0 is then an inflation of V
Params.acap=0;
Params.zcap=(z_grid(end-1)+z_grid(end))/2; % between the top two z: the cap bites at z_top only

a_grid=5*linspace(0,1,101)'.^3; % the bank's usual grid; a_grid(1)=0=acap is the only alive node at z_top
n_a=length(a_grid);

pi_z(1:end-2,end)=0; % the three lowest z never reach z_top
pi_z=pi_z./sum(pi_z,2);
fprintf('\ncap at top z, semiz: V(a,semiz,z_top)=-Inf for a>%g; pi_z column %i zeroed in rows 1:%i and rows renormalised \n',Params.acap,n_z,n_z-2);

vfoptions=struct();
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
vfoptions_e=vfoptions;
vfoptions_e.n_e=vfoptionsbaseline.n_e;
vfoptions_e.e_grid=vfoptionsbaseline.e_grid;
vfoptions_e.pi_e=vfoptionsbaseline.pi_e;

%% nod1, z, noe
ReturnFn=@(d2,aprime,a,semiz,z,r,w,sigma,uempbenefit,searcheffortcost,acap,zcap) ReturnFn_nod1_z_noe_InfNaNcap_semiz(d2,aprime,a,semiz,z,r,w,sigma,uempbenefit,searcheffortcost,acap,zcap);
output=CoreInfHorzVFIAlgo_InfNaNcompare('cap at top z: nod1, z, noe, semiz',n_d2_semiz,n_a,n_z,d2_grid_semiz,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions,figure_c);

%% d1, z, noe
ReturnFn=@(d1,d2,aprime,a,semiz,z,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost,acap,zcap) ReturnFn_d1_z_noe_InfNaNcap_semiz(d1,d2,aprime,a,semiz,z,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost,acap,zcap);
output=CoreInfHorzVFIAlgo_InfNaNcompare('cap at top z: d1, z, noe, semiz',n_d_semiz,n_a,n_z,d_grid_semiz,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions,figure_c+1);

%% nod1, z, e
ReturnFn=@(d2,aprime,a,semiz,z,e,r,w,sigma,uempbenefit,searcheffortcost,acap,zcap) ReturnFn_nod1_z_e_InfNaNcap_semiz(d2,aprime,a,semiz,z,e,r,w,sigma,uempbenefit,searcheffortcost,acap,zcap);
output=CoreInfHorzVFIAlgo_InfNaNcompare('cap at top z: nod1, z, e, semiz',n_d2_semiz,n_a,n_z,d2_grid_semiz,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions_e,figure_c+2);

%% d1, z, e
ReturnFn=@(d1,d2,aprime,a,semiz,z,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost,acap,zcap) ReturnFn_d1_z_e_InfNaNcap_semiz(d1,d2,aprime,a,semiz,z,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost,acap,zcap);
output=CoreInfHorzVFIAlgo_InfNaNcompare('cap at top z: d1, z, e, semiz',n_d_semiz,n_a,n_z,d_grid_semiz,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions_e,figure_c+3);

end
