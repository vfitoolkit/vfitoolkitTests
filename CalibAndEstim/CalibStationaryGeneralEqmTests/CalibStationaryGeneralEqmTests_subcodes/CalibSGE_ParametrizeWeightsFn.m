function Parameters=CalibSGE_ParametrizeWeightsFn(Parameters)
% ParametrizeParamsFn used by CalibStationaryGeneralEqmTests parts 15 and 20: the type masses (PTypeDistParamNames, ptypemass) from one
% calibrated scalar, the share of the first type, so that the masses always sum to one while the share is being calibrated.

Parameters.ptypemass=[Parameters.shareA; 1-Parameters.shareA];

end
