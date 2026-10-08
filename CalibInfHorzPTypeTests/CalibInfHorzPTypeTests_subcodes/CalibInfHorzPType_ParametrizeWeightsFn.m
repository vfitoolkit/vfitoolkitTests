function Parameters=CalibInfHorzPType_ParametrizeWeightsFn(Parameters)
% ParametrizeParamsFn used by CalibInfHorzPTypeTests part 4 (l): the type weights (PTypeDistParamNames, ptypemass) from one calibrated
% scalar, the share of the low type, so that the weights always sum to one while the share is being calibrated.

Parameters.ptypemass=[Parameters.sharelow; 1-Parameters.sharelow];

end
