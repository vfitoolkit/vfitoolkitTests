function [z_grid,pi_z]=CalibFHorz_ExogShockFn(rho_z,sigma_z)
% ExogShockFn used by CalibFHorzTests part 5: the bank's 5-state AR(1) in logs, by Farmer-Toda, with the autocorrelation and
% innovation standard deviation as parameters (so rho_z can be calibrated). n_z=5 is hardcoded as ExogShockFn inputs can only be
% parameters; at rho_z=0.9, sigma_z=0.03 this is exactly the z_grid and pi_z of the bank's driver.

[z_grid,pi_z]=discretizeAR1_FarmerToda(0,rho_z,sigma_z,5);
z_grid=exp(z_grid);

end
