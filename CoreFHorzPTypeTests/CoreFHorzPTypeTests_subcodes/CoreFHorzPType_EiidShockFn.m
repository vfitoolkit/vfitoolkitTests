function [e_grid,pi_e]=CoreFHorzPType_EiidShockFn(sigma_e)
% EiidShockFn used by CoreFHorzPType_ShockFns: a 3-state iid shock in logs, by Farmer-Toda.
% n_e=3 is hardcoded as EiidShockFn inputs can only be parameters.

[e_grid,pi_e]=discretizeAR1_FarmerToda(0,0,sigma_e,3);
pi_e=pi_e(1,:)';
e_grid=exp(e_grid);

end
