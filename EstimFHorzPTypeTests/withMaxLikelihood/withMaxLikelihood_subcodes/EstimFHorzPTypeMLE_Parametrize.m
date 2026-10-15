function Params=EstimFHorzPTypeMLE_Parametrize(Params)
% ParametrizeParamsFn of the PType maximum likelihood tests: the masses of the two permanent types from the share of the second
% (so that estimating pshare_high keeps the masses summing to one)
Params.ptypeweights=[1-Params.pshare_high; Params.pshare_high];
end
