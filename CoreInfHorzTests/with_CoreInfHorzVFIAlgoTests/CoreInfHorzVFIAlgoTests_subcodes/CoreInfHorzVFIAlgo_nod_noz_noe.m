function output=CoreInfHorzVFIAlgo_nod_noz_noe(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)

vfoptions=struct();
n_d=0; d_grid=[];             % without d
n_z=0; z_grid=[]; pi_z=[];    % without z
% without e (do not set vfoptions.n_e)

ReturnFn=@(aprime,a,r,w,sigma) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,sigma);

output=CoreInfHorzVFIAlgo_algocompare('nod, noz, noe',n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);

end
