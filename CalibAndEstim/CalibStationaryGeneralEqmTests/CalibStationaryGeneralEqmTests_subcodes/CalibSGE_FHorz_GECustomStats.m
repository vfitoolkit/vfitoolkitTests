function CustomStats=CalibSGE_FHorz_GECustomStats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,heteroagentoptions,vfoptions,simoptions)
% heteroagentoptions.CustomModelStats used by CalibStationaryGeneralEqmTests part 5: Kcustom, aggregate capital computed by hand from
% the agent distribution, which the CapitalMarket condition then uses in place of the aggregate K (so the GE conditions are the
% baseline ones). V is part of the interface (and optionsflag must be 0: HeteroAgentStationaryEqm_InfHorz used to hand CustomModelStats its internal vfoptions): the joint objective only kept it when caliboptions.CustomModelStats was set.

N_a=prod(n_a);
N_z=prod(n_z);
StationaryDist=reshape(StationaryDist,[N_a,N_z*N_j]);
CustomStats.Kcustom=sum(sum(StationaryDist.*a_grid));
CustomStats.optionsflag=double(isfield(vfoptions,'alreadygridvals') && vfoptions.alreadygridvals==1); % 1 if handed the internal vfoptions; it is added to CapitalMarket, so the GE conditions only match the baseline when the user's own vfoptions are handed in

end
