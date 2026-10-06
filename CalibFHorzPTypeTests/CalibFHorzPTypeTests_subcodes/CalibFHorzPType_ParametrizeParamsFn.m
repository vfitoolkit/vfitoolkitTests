function Parameters=CalibFHorzPType_ParametrizeParamsFn(Parameters)
% ParametrizeParamsFn used by CalibFHorzPTypeTests part 4: the per-type earnings profile kappa_j_pt (N_i by N_j) is a straight line
% in age with a per-type intercept and a common slope, kappa_j_pt(ii,:) = kappa0.(type) + kappaslope*(agej-1) while working and zero
% in retirement, so that the structure kappa0 (one value per type) and the scalar kappaslope can be calibrated in place of kappa_j_pt.
% The calibration calls this on every evaluation, after setting the calibrated parameters, and before solving the model.

Parameters.kappa_j_pt=[Parameters.kappa0.low+Parameters.kappaslope*(Parameters.agej-1); Parameters.kappa0.high+Parameters.kappaslope*(Parameters.agej-1)].*(Parameters.agej<Parameters.Jr);

end
