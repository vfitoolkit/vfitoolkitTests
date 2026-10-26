function CustomStats=CalibFHorz_CustomModelStats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions)
% CustomModelStats used by CalibFHorzTests part 5. Two statistics:
%   meanassets: the mean of assets computed by hand from the agent distribution (so it does not depend on which form z_grid arrives in);
%               its target is the AllStats mean of assets at the truth.
%   zgridcheck: max|z_grid - caliboptions.CalibFHorz_expected_zgrid| (1 if the sizes differ); its target is 0. With
%               caliboptions.CustomModelStats_usergrids=1 the expected grid is the user's z_grid, with =0 the internal joint grid
%               z_gridvals_J, so this checks which grid the calibration hands to CustomModelStats.
% (V, Policy, FnsToEvaluate, n_d, d_grid, pi_z, vfoptions, simoptions are part of the interface and unused here.)

N_a=prod(n_a);
N_z=prod(n_z);
StationaryDist=reshape(StationaryDist,[N_a,N_z*N_j]);
CustomStats.meanassets=gather(sum(sum(StationaryDist.*a_grid)));

expected=caliboptions.CalibFHorz_expected_zgrid;
if numel(z_grid)==numel(expected)
    CustomStats.zgridcheck=gather(max(abs(z_grid(:)-expected(:))));
else
    CustomStats.zgridcheck=1;
end

end
