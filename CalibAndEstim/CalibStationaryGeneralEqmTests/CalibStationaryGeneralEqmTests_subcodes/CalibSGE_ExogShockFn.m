function [z_grid,pi_z]=CalibSGE_ExogShockFn(rho_z,sigma_z)
% ExogShockFn used by CalibStationaryGeneralEqmTests parts 5 and 11: the bank's 5-state AR(1) in logs, by Farmer-Toda, normalised to
% E[z]=1, with the autocorrelation and innovation standard deviation as parameters (so rho_z can be calibrated). n_z=5 is hardcoded
% as ExogShockFn inputs can only be parameters; at rho_z=0.9 and sigma_z=0.03 (FHorz) or 0.2 (InfHorz) this is exactly the z_grid and pi_z of the bank's driver.

[z_grid,pi_z]=discretizeAR1_FarmerToda(0,rho_z,sigma_z,5);
z_grid=exp(z_grid);
pi_z_stat=pi_z^1000; pi_z_stat=pi_z_stat(1,:)';
z_grid=z_grid./sum(z_grid.*pi_z_stat);

end
