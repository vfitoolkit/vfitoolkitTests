function Params=EstimParams_to_Params(Params,EstimParams)
% EstimFHorzTests: copies the estimated parameters (a structure, as the estimation commands return them) into the parameter structure.
% With permanent types a per-type estimate is a structure with a field per type, which is copied as it is.
names=fieldnames(EstimParams);
for pp=1:length(names)
    if ~strcmp(names{pp},'notes')
        Params.(names{pp})=EstimParams.(names{pp});
    end
end
end
