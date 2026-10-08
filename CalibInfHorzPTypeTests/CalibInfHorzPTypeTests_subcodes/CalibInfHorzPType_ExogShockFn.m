function [z_grid,pi_z]=CalibInfHorzPType_ExogShockFn(rho_z,sigma_z)
% ExogShockFn used by CalibInfHorzPTypeTests part 5, common or per type (then wrapped as vfoptions.ExogShockFn.low=@(rho_z_low,sigma_z) ...):
% the bank's 5-state AR(1) in logs, by Farmer-Toda. n_z=5 is hardcoded as ExogShockFn inputs can only be parameters; at rho_z=0.9,
% sigma_z=0.1 this is exactly the z_grid and pi_z of the bank's driver.

[z_grid,pi_z]=discretizeAR1_FarmerToda(0,rho_z,sigma_z,5);
z_grid=exp(z_grid);

end
