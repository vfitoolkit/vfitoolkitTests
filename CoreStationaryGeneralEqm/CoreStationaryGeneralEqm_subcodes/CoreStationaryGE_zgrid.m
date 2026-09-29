function [z_grid,pi_z]=CoreStationaryGE_zgrid(rho_z,sigma_epsz,n_z,r,rgearing)
% The z process of CoreStationaryGE_setup, written as a function of parameters so that it can
% be handed to the toolkit as vfoptions.ExogShockFn. Because one of its inputs is a general eqm
% price, doing so makes the shock grids a thing the general eqm solver has to rebuild at every
% price vector it tries (it sets heteroagentoptions.gridsinGE=1 internally).
%
% rgearing is what separates the two variants of that test:
%   rgearing=0  r is an input, so the toolkit treats the grids as general-eqm-determined and
%               rebuilds them on every iteration, but what it rebuilds is EXACTLY the baseline
%               grid. The model is therefore unchanged, and the general eqm conditions must come
%               out bit-identical to the baseline run: any difference at all is the rebuild
%               machinery corrupting something rather than economics.
%   rgearing>0  income risk genuinely rises with r, so the grids really do move with the price.
%               The equilibrium then changes, and what gets checked instead is that the grids the
%               solver was using at its own answer are the grids this function gives at that
%               answer - i.e. that the rebuild tracks the current prices instead of being frozen
%               at the initial guess. See the frozen-grid checks in the extraoptions subcodes.
%
% This must reproduce CoreStationaryGE_setup line for line when rgearing=0, normalization of
% E[z] to one included, or the neutral check is comparing two different models and proves nothing.
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,rho_z,sigma_epsz+rgearing*r,n_z);
z_grid=exp(z_grid);
% normalize so E[z]=1 under the stationary distribution
pi_z_stat=pi_z^1000; pi_z_stat=pi_z_stat(1,:)';
z_grid=z_grid./sum(z_grid.*pi_z_stat);

end
