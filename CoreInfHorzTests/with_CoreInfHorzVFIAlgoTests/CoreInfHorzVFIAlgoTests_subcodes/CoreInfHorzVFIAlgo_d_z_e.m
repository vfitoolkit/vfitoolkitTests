function output=CoreInfHorzVFIAlgo_d_z_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)

vfoptions=struct();
% with d
% with z (markov)
% with e (iid)
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;

ReturnFn=@(d,aprime,a,z,e,r,w,sigma,eta,varphi) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,sigma,eta,varphi);

output=CoreInfHorzVFIAlgo_algocompare('d, z, e',n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);

end
