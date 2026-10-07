function [z_grid,pi_z]=CoreInfHorzPType_ExogShockFn(rho_z,sigma_z)
% The exogenous shock of CoreInfHorzPType_ShockFns, as a function of its parameters: a 5-state AR(1) in logs (Farmer-Toda).
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,rho_z,sigma_z,5);
z_grid=exp(z_grid);
end
