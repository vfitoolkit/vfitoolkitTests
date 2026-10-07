function [e_grid,pi_e]=CalibFHorz_EiidShockFn(sigma_e)
% EiidShockFn used by CalibFHorzTests part 7: a 3-state iid shock in logs, by Farmer-Toda, with its standard deviation as a
% parameter (so sigma_e can be calibrated). n_e=3 is hardcoded as EiidShockFn inputs can only be parameters.

[e_grid,pi_e]=discretizeAR1_FarmerToda(0,0,sigma_e,3);
pi_e=pi_e(1,:)';
e_grid=exp(e_grid);

end
