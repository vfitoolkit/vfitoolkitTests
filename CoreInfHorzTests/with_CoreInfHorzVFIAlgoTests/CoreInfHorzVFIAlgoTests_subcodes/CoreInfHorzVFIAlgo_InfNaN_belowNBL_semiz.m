function output=CoreInfHorzVFIAlgo_InfNaN_belowNBL_semiz(n_d2_semiz,n_d_semiz,d2_grid_semiz,d_grid_semiz,n_z,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,figure_c)
% Deliberate trigger, device 1 (the asset grid extends below the natural borrowing limit), for the models with a
% semi-exogenous state semiz (employment, its transitions depend on the search effort d2). See
% CoreInfHorzVFIAlgo_InfNaN_belowNBL for the device; the dead region is the same (nodes below the limit at every
% (semiz,z), since an alive a' has to be alive at every reachable (semiz',z')) and the unemployed (income uempbenefit=0.2,
% limit -4 at r=0.05) are the ones whose constraint binds first.
%
% Runs the eight semiz variants: nod1/d1, z/noz, e/noe (the semiz e raws do NOT have the per-term form yet, unlike the
% plain e raws, so here e is not a positive control). The SemiExo postGI raws, their EndogenousVFI twins (the endoVFI
% configs) and, without z, their RelativeVFI twins. Eight figures, figure_c to figure_c+7.

DF=DiscountFactorParamNames;
Params.w=0.5; % so that c<1, u<0, at every state that matters: a continuation priced at 0 is then an inflation of V

n_neg=12;
a_grid=[(-n_neg:1:-1)'; 5*linspace(0,1,101)'.^3];
n_a=length(a_grid);
fprintf('\nbelow NBL, semiz: a_grid has %i integer-spaced nodes from %i to -1 below the usual cubic grid on [0,5] (n_a=%i); natural borrowing limit when unemployed, uempbenefit=%g, r=%g: %.2f \n',n_neg,-n_neg,n_a,Params.uempbenefit,Params.r,-Params.uempbenefit/Params.r);

vfoptions=struct();
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
vfoptions_e=vfoptions;
vfoptions_e.n_e=vfoptionsbaseline.n_e;
vfoptions_e.e_grid=vfoptionsbaseline.e_grid;
vfoptions_e.pi_e=vfoptionsbaseline.pi_e;

%% nod1, z, noe (SemiExo_postGI_nod1_raw, SemiExo_postGI_sparse_nod1_raw, EndogenousVFI_SemiExo_postGI_nod1_raw)
ReturnFn=@(d2,aprime,a,semiz,z,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,sigma,uempbenefit,searcheffortcost);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: nod1, z, noe, semiz',n_d2_semiz,n_a,n_z,d2_grid_semiz,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions,figure_c);

%% d1, z, noe (SemiExo_Refine_postGI_raw, EndogenousVFI_SemiExo_Refine_postGI_raw)
ReturnFn=@(d1,d2,aprime,a,semiz,z,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_z_noe_semiz(d1,d2,aprime,a,semiz,z,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: d1, z, noe, semiz',n_d_semiz,n_a,n_z,d_grid_semiz,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions,figure_c+1);

%% nod1, z, e (SemiExo_postGI_nod1_e_raw, SemiExo_postGI_sparse_nod1_e_raw, EndogenousVFI_SemiExo_postGI_nod1_e_raw)
ReturnFn=@(d2,aprime,a,semiz,z,e,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_z_e_semiz(d2,aprime,a,semiz,z,e,r,w,sigma,uempbenefit,searcheffortcost);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: nod1, z, e, semiz',n_d2_semiz,n_a,n_z,d2_grid_semiz,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions_e,figure_c+2);

%% d1, z, e (SemiExo_Refine_postGI_e_raw, EndogenousVFI_SemiExo_Refine_postGI_e_raw)
ReturnFn=@(d1,d2,aprime,a,semiz,z,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_z_e_semiz(d1,d2,aprime,a,semiz,z,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: d1, z, e, semiz',n_d_semiz,n_a,n_z,d_grid_semiz,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions_e,figure_c+3);

%% nod1, noz, noe (SemiExo_postGI_nod1_noz_raw, sparse twin, RelativeVFI_SemiExo_postGI_nod1_noz_raw)
ReturnFn=@(d2,aprime,a,semiz,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_semiz(d2,aprime,a,semiz,r,w,sigma,uempbenefit,searcheffortcost);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: nod1, noz, noe, semiz',n_d2_semiz,n_a,0,d2_grid_semiz,a_grid,[],[],ReturnFn,Params,DF,vfoptions,figure_c+4);

%% d1, noz, noe (SemiExo_Refine_postGI_noz_raw, RelativeVFI_SemiExo_Refine_postGI_noz_raw)
ReturnFn=@(d1,d2,aprime,a,semiz,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_noz_noe_semiz(d1,d2,aprime,a,semiz,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: d1, noz, noe, semiz',n_d_semiz,n_a,0,d_grid_semiz,a_grid,[],[],ReturnFn,Params,DF,vfoptions,figure_c+5);

%% nod1, noz, e (SemiExo_postGI_nod1_noz_e_raw, sparse twin, RelativeVFI_SemiExo_postGI_nod1_noz_e_raw)
ReturnFn=@(d2,aprime,a,semiz,e,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_e_semiz(d2,aprime,a,semiz,e,r,w,sigma,uempbenefit,searcheffortcost);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: nod1, noz, e, semiz',n_d2_semiz,n_a,0,d2_grid_semiz,a_grid,[],[],ReturnFn,Params,DF,vfoptions_e,figure_c+6);

%% d1, noz, e (SemiExo_Refine_postGI_noz_e_raw, RelativeVFI_SemiExo_Refine_postGI_noz_e_raw)
ReturnFn=@(d1,d2,aprime,a,semiz,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_noz_e_semiz(d1,d2,aprime,a,semiz,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
output=CoreInfHorzVFIAlgo_InfNaNcompare('below NBL: d1, noz, e, semiz',n_d_semiz,n_a,0,d_grid_semiz,a_grid,[],[],ReturnFn,Params,DF,vfoptions_e,figure_c+7);

end
