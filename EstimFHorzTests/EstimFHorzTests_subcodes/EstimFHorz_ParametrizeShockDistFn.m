function Parameters=EstimFHorz_ParametrizeShockDistFn(Parameters)
% ParametrizeParamsFn used by EstimFHorzTests part 9: the estimated primitives are half_rho and loginitdecay; the DERIVED parameters
% rho_z (which feeds vfoptions.ExogShockFn) and initdecay (which feeds the jequaloneDist function) are written here. Neither
% estimated primitive appears in any shock-fn or dist-fn signature, so the estimation command can only detect the dependence through
% ParseParametrizeParamsFn (or a user-given estimoptions.ParametrizeParamsFn_OutputNames) -- which is what part 9 pins down (before
% the 2026-10-12 fix the shocks were silently frozen at the initial parameters in this situation).

Parameters.rho_z=2*Parameters.half_rho;
Parameters.initdecay=exp(Parameters.loginitdecay);

end
