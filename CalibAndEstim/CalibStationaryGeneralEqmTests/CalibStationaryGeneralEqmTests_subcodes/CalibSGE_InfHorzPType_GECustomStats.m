function CustomStats=CalibSGE_InfHorzPType_GECustomStats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,heteroagentoptions,vfoptions,simoptions)
% heteroagentoptions.CustomModelStats used by CalibStationaryGeneralEqmTests part 21 (CalibrateBIHAModel_PType, and the inner solve of
% HeteroAgentStationaryEqm_InfHorz_PType): Kcustom, aggregate capital computed by hand from the agent distribution (each type's weighted
% by its mass), which the CapitalMarket condition then uses in place of the aggregate K (so the GE conditions are the baseline ones).
% optionsflag is 1 if handed the internal vfoptions (alreadygridvals=1); it is added to CapitalMarket, so the GE conditions only match
% the baseline when the user's own vfoptions are handed in. V is part of the interface.

N_a=prod(n_a);
Kcustom=0;
for ii=1:length(Names_i)
    dist_ii=StationaryDist.(Names_i{ii});
    dist_ii=reshape(dist_ii,[N_a,numel(dist_ii)/N_a]);
    Kcustom=Kcustom+StationaryDist.ptweights(ii)*sum(sum(dist_ii.*a_grid));
end
CustomStats.Kcustom=Kcustom;
CustomStats.optionsflag=double(isfield(vfoptions,'alreadygridvals') && vfoptions.alreadygridvals==1);

end
