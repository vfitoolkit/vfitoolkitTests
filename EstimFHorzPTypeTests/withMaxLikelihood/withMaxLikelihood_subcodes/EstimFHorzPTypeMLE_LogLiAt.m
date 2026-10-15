function [LogLi,Post]=EstimFHorzPTypeMLE_LogLiAt(vec,Params,Data,MeasurementError,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,DiscountFactorParamNames,jequaloneDist,PTypeDistParamNames,vfoptions,simoptions)
% EstimFHorzPTypeMLETests doPart(2): a direct solve and log-likelihood at the parameter vector of the estimation, in its block order
% [beta; sig_c of the first type; sig_c of the second type; pshare_high], with the masses of the types from EstimFHorzPTypeMLE_Parametrize.
Params.beta=vec(1);
Params.sig_c.(Names_i{1})=vec(2);
Params.sig_c.(Names_i{2})=vec(3);
Params.pshare_high=vec(4);
Params=EstimFHorzPTypeMLE_Parametrize(Params);
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
[~,LogLi,Post]=LogLikelihood_FHorz_PType(Data,MeasurementError,jequaloneDist,PTypeDistParamNames,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions);
LogLi=gather(LogLi);
end
