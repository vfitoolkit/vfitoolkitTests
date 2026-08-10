function output=CoreInfHorzVFIAlgo_d_noz_noe(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)

vfoptions=struct();
% with d
n_z=0; z_grid=[]; pi_z=[];    % without z
% without e (do not set vfoptions.n_e)

ReturnFn=@(d,aprime,a,r,w,sigma,eta,varphi) ReturnFn_d_noz_noe_nosemiz(d,aprime,a,r,w,sigma,eta,varphi);

output=CoreInfHorzVFIAlgo_algocompare('d, noz, noe',n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);

end
