function Parameters=CalibInfHorzPType_ParametrizeParamsFn(Parameters)
% ParametrizeParamsFn used by CalibInfHorzPTypeTests part 4: a common vector parameter theta=[log w_low; log w_high; discount rate] is
% mapped to the per-type wage w_pt (as a structure over the types) and the discount factor beta=1/(1+theta(3)), so that one vector
% can be calibrated, in whole or in part through caliboptions.omitcalibparam, in place of the per-type wage and beta. The calibration
% calls this on every evaluation, after setting the calibrated parameters, and before solving the model.

Parameters.w_pt=struct();
Parameters.w_pt.low=exp(Parameters.theta(1));
Parameters.w_pt.high=exp(Parameters.theta(2));
Parameters.beta=1/(1+Parameters.theta(3));

end
