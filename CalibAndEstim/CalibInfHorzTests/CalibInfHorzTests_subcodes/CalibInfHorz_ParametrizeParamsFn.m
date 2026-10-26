function Parameters=CalibInfHorz_ParametrizeParamsFn(Parameters)
% ParametrizeParamsFn used by CalibInfHorzTests part 4: a vector parameter theta=[log w; discount rate] is mapped to the wage
% w=exp(theta(1)) and the discount factor beta=1/(1+theta(2)), so that a vector parameter can be calibrated (in whole, or in part
% through caliboptions.omitcalibparam) in place of the two scalars. The calibration calls this on every evaluation, after setting
% the calibrated parameters, and before solving the model.

Parameters.w=exp(Parameters.theta(1));
Parameters.beta=1/(1+Parameters.theta(2));

end
