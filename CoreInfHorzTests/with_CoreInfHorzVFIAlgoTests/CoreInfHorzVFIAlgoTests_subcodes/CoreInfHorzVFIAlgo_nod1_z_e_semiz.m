function output=CoreInfHorzVFIAlgo_nod1_z_e_semiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)
% With a semi-exogenous state semiz (its transitions depend on the decision variable d2), without d1, with z, with an iid e.

vfoptions=struct();
% with semiz
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
% with e (iid)
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;

ReturnFn=@(d2,aprime,a,semiz,z,e,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_z_e_semiz(d2,aprime,a,semiz,z,e,r,w,sigma,uempbenefit,searcheffortcost);

output=CoreInfHorzVFIAlgo_algocompare('nod1, z, e, semiz',n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);

end
