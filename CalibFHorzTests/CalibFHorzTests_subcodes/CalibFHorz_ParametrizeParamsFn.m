function Parameters=CalibFHorz_ParametrizeParamsFn(Parameters)
% ParametrizeParamsFn used by CalibFHorzTests part 4: the deterministic earnings profile kappa_j is a straight line in age,
% kappa_j = kappa0 + kappaslope*(agej-1) while working and zero in retirement, so that the two scalars kappa0 and kappaslope can be
% calibrated in place of the age-dependent kappa_j. The calibration calls this on every evaluation, after setting the calibrated
% parameters, and before solving the model.

Parameters.kappa_j=(Parameters.kappa0+Parameters.kappaslope*(Parameters.agej-1)).*(Parameters.agej<Parameters.Jr);

end
