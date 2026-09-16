function output=DiscP4_downstream(calib)
% P4: feed this block's outputs into an actual toolkit solve
%
% This block matters more than the other downstream subcodes, because it is the only one whose
% output is a STACKED two-dimensional grid with a joint transition matrix. That is a shape the
% solvers have to be told about, via n_z=[xnum,znum], and getting it wrong is exactly the kind of
% error that is invisible in the discretization tests and obvious the moment something consumes it.

fprintf('\n========== P4: downstream (do the outputs go into a solve) ========== \n')

output=struct();
xnum=5; znum=7;
n_d=0; d_grid=[];
n_a=5; a_grid=linspace(0.1,5,n_a)';
Params=struct(); Params.beta=0.9; Params.r=0.03;
DiscountFactorParamNames={'beta'};
% Two exogenous variables now, so the ReturnFn takes both. Only z enters the budget; x is the log
% volatility and affects the model only through the transition matrix.
ReturnFn=@(aprime,a,x,z,r) log(max((1+r)*a+exp(z)-aprime,10^(-10)));
FnsToEvaluate.zvalue=@(aprime,a,x,z) z;
FnsToEvaluate.xvalue=@(aprime,a,x,z) x;

cmdname={'AR1wSV_FarmerToda','AR1wSV_Tauchen'};
for c_c=1:2
    if c_c==1
        opts=struct(); opts.nSigmas=3;
        [z_grid,pi_z]=discretizeAR1wSV_FarmerToda(calib.rho,calib.phi,calib.sigmau,calib.sigmae,xnum,znum,opts);
    else
        [z_grid,pi_z]=discretizeAR1wSV_Tauchen(calib.rho,calib.phi,calib.sigmau,calib.sigmae,xnum,znum,3,struct());
    end
    n_z=[xnum,znum]; % the stacked convention: x first, x fastest

    vfoptions=struct(); simoptions=struct();
    [V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params,[]);
    AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

    fprintf('%s: solve runs, V has the right size [T0], this should be zero: %i \n',cmdname{c_c},any(size(V)~=[n_a,xnum,znum])&&any(size(V)~=[n_a,xnum*znum]))
    fprintf('%s: StationaryDist sums to one [T1], this should be zero: %2.8e \n',cmdname{c_c},abs(sum(StationaryDist(:))-1))
    fprintf('%s: V is finite everywhere [T0], this should be zero: %i \n',cmdname{c_c},any(~isfinite(V(:))))
    % The substantive one: the toolkit's own mean of x must match the chain's, which it can only do
    % if the stacked grid was unpacked the way the discretization packed it.
    fprintf('%s: AllStats mean of x is %2.4f, against xBar %2.4f [T2], this should be small: %2.3e \n',cmdname{c_c},AllStats.xvalue.Mean,calib.xBar,abs(AllStats.xvalue.Mean-calib.xBar))
    fprintf('%s: AllStats mean of z is %2.3e, against zero [T2], this should be small: %2.3e \n',cmdname{c_c},AllStats.zvalue.Mean,abs(AllStats.zvalue.Mean))
    output.cmd(c_c).name=cmdname{c_c};
end

end
