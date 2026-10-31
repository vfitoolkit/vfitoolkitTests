function CustomStats=CalibSGE_InfHorzPType_CustomModelStats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions)
% caliboptions.CustomModelStats used by CalibStationaryGeneralEqmTests part 21 (CalibrateBIHAModel_PType). Three statistics:
%   meanK:      the mean of assets computed by hand from the agent distribution (each type's distribution weighted by its mass,
%               StationaryDist.ptweights); its target is the grouped AllStats mean of assets at the truth.
%   zgridcheck: max|z_grid - caliboptions.CalibSGE_expected_zgrid| for the first type (z_grid comes as a struct over the types; 1 if the
%               sizes differ); its target is 0, so this checks which grid (the user's, or the internal one) the calibration hands in.
%   optionscheck: 1 if the vfoptions handed in have alreadygridvals=1 (the internal options), 0 otherwise; its target is 0.
% (V, Policy, FnsToEvaluate, n_d, d_grid, pi_z, simoptions are part of the interface and unused here.)

N_a=prod(n_a);
meanK=0;
for ii=1:length(Names_i)
    dist_ii=StationaryDist.(Names_i{ii});
    dist_ii=reshape(dist_ii,[N_a,numel(dist_ii)/N_a]);
    meanK=meanK+StationaryDist.ptweights(ii)*sum(sum(dist_ii.*a_grid));
end
CustomStats.meanK=gather(meanK);

if isstruct(z_grid)
    zg=z_grid.(Names_i{1});
else
    zg=z_grid;
end
expected=caliboptions.CalibSGE_expected_zgrid;
if numel(zg)==numel(expected)
    CustomStats.zgridcheck=gather(max(abs(zg(:)-expected(:))));
else
    CustomStats.zgridcheck=1;
end
CustomStats.optionscheck=double(isfield(vfoptions,'alreadygridvals') && vfoptions.alreadygridvals==1);

end
