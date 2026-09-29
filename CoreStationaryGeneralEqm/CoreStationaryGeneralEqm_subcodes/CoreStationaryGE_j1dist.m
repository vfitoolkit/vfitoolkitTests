function jequaloneDist=CoreStationaryGE_j1dist(a_grid,z_grid,n_a,n_z)
% jequaloneDist as a FUNCTION rather than as a matrix: newborns all hold the lowest asset grid
% point, and their productivity weights are proportional to the productivity grid itself.
%
% Why the weights are read off z_grid, rather than being the stationary distribution of pi_z that
% CoreStationaryGE_setup uses for the matrix form: a jequaloneDist function is handed
% (a_grid,z_grid,n_a,n_z) and NEVER pi_z, so a newborn distribution defined from the transition
% matrix cannot be written in this form at all. Defining it from z_grid is also what makes the test
% sharp - it comes out wrong if the toolkit hands over the wrong z_grid, either the internal
% joint-grid form instead of the user's own, or a grid that was not rebuilt at the current prices.
%
% a_grid is part of that same fixed interface, and is deliberately unused here.
jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,:)=shiftdim(z_grid./sum(z_grid),-1);

end
