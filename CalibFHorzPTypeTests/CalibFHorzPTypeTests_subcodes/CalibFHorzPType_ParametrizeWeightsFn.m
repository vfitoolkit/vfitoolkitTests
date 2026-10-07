function Parameters=CalibFHorzPType_ParametrizeWeightsFn(Parameters)
% ParametrizeParamsFn used by CalibFHorzPTypeTests part 7: the type weights (PTypeDistParamNames) from one calibrated scalar, the
% share of the low type, so that the weights always sum to one while the share is being calibrated.

Parameters.ptypeweights=[Parameters.sharelow; 1-Parameters.sharelow];

end
