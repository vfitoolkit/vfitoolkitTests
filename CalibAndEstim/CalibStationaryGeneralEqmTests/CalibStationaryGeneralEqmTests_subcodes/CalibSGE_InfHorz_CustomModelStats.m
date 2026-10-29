function CustomStats=CalibSGE_InfHorz_CustomModelStats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions)
% caliboptions.CustomModelStats used by CalibStationaryGeneralEqmTests part 11. Three statistics:
%   meanK:      the mean of assets computed by hand from the agent distribution; its target is the AllStats mean of assets at the truth.
%   zgridcheck: max|z_grid - caliboptions.CalibSGE_expected_zgrid| (1 if the sizes differ); its target is 0, so this checks which grid
%               (the user's, or the internal one) the calibration hands to CustomModelStats.
%   optionscheck: 1 if the vfoptions handed in have alreadygridvals=1 (the internal options), 0 otherwise; its target is 0 (the
%               calibration must hand CustomModelStats the user's own vfoptions, as CalibrateOLGModel does).
% (V, Policy, FnsToEvaluate, n_d, d_grid, pi_z, simoptions are part of the interface and unused here.)

N_a=prod(n_a);
StationaryDist=reshape(StationaryDist,[N_a,numel(StationaryDist)/N_a]);
CustomStats.meanK=gather(sum(sum(StationaryDist.*a_grid)));

expected=caliboptions.CalibSGE_expected_zgrid;
if numel(z_grid)==numel(expected)
    CustomStats.zgridcheck=gather(max(abs(z_grid(:)-expected(:))));
else
    CustomStats.zgridcheck=1;
end
CustomStats.optionscheck=double(isfield(vfoptions,'alreadygridvals') && vfoptions.alreadygridvals==1);

end
