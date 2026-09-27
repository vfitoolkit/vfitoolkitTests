function [z_grid,pi_z]=CoreFHorzPType_ExogShockFn(rho_z,sigma_z)
% ExogShockFn used by CoreFHorzPType_ShockFns: a 5-state AR(1) in logs, by Farmer-Toda.
% n_z=5 is hardcoded as ExogShockFn inputs can only be parameters.

[z_grid,pi_z]=discretizeAR1_FarmerToda(0,rho_z,sigma_z,5);
z_grid=exp(z_grid);

end
