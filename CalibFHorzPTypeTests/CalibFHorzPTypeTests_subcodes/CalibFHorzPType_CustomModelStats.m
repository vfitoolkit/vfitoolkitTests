function CustomStats=CalibFHorzPType_CustomModelStats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions)
% CustomModelStats used by CalibFHorzPTypeTests part 5 (the PType signature has Names_i). Two statistics:
%   meanassets: the mean of assets over the whole population, by hand from the per-type agent distributions and the type weights
%               (so it does not depend on which form z_grid arrives in); its target is the grouped AllStats mean of assets at the truth.
%   zgridcheck: max|z_grid - caliboptions.CalibFHorz_expected_zgrid| (1 if the shapes differ), over the fields when the grids are a
%               structure over the types; its target is 0. With caliboptions.CustomModelStats_usergrids=1 the expected grid is the
%               user's z_grid, with =0 the internal form that ExogShockSetup_FHorz_PType returns.
% (V, Policy, FnsToEvaluate, n_d, N_j, d_grid, pi_z, vfoptions, simoptions are part of the interface and unused here.)

N_a=prod(n_a);
N_z=prod(n_z);
meanassets=0;
for ii=1:length(Names_i)
    dist_ii=reshape(StationaryDist.(Names_i{ii}),[N_a,numel(StationaryDist.(Names_i{ii}))/N_a]);
    meanassets=meanassets+StationaryDist.ptweights(ii)*sum(sum(dist_ii.*a_grid));
end
CustomStats.meanassets=gather(meanassets);

expected=caliboptions.CalibFHorz_expected_zgrid;
if isstruct(z_grid) && isstruct(expected)
    d=0;
    fn=fieldnames(expected);
    for ii=1:length(fn)
        if isfield(z_grid,fn{ii}) && numel(z_grid.(fn{ii}))==numel(expected.(fn{ii}))
            d=max(d,gather(max(abs(z_grid.(fn{ii})(:)-expected.(fn{ii})(:)))));
        else
            d=1;
        end
    end
    CustomStats.zgridcheck=d;
elseif ~isstruct(z_grid) && ~isstruct(expected) && numel(z_grid)==numel(expected)
    CustomStats.zgridcheck=gather(max(abs(z_grid(:)-expected(:))));
else
    CustomStats.zgridcheck=1;
end

end
