function output=CoreInfHorzVFIAlgo_nod_noz_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)

vfoptions=struct();
n_d=0; d_grid=[];             % without d
n_z=0; z_grid=[]; pi_z=[];    % without z
% with e (iid)
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;

ReturnFn=@(aprime,a,e,r,w,sigma) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,sigma);

output=CoreInfHorzVFIAlgo_algocompare('nod, noz, e',n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);

end
