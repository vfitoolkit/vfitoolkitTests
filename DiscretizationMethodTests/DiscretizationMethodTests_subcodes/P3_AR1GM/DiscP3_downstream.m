function output=DiscP3_downstream(calib)
% P3: feed this block's outputs into an actual toolkit solve
%
% Same minimal purpose as the other blocks': the solve runs, the shapes are right, and the z
% marginal of the agent distribution is the chain's own stationary distribution.

fprintf('\n========== P3: downstream (do the outputs go into a solve) ========== \n')

output=struct();
znum=9;
n_d=0; d_grid=[];
n_a=5; a_grid=linspace(0.1,5,n_a)';
Params=struct(); Params.beta=0.9; Params.r=0.03;
DiscountFactorParamNames={'beta'};
ReturnFn=@(aprime,a,z,r) log(max((1+r)*a+exp(z)-aprime,10^(-10)));
FnsToEvaluate.zvalue=@(aprime,a,z) z;

cmdname={'AR1wGM_FarmerToda','AR1wGM_Tauchen'};
for c_c=1:2
    if c_c==1
        [z_grid,pi_z,~]=discretizeAR1wGM_FarmerToda(calib.mew,calib.rho,calib.mixprobs,calib.mu,calib.sigma,znum,struct('verbose',0));
    else
        [z_grid,pi_z]=discretizeAR1wGM_Tauchen(calib.mew,calib.rho,calib.mixprobs,calib.mu,calib.sigma,znum,3,struct());
    end
    n_z=znum;
    vfoptions=struct(); simoptions=struct();
    [V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params,[]);
    AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

    fprintf('%s: solve runs, size of V [T0], this should be zero: %i \n',cmdname{c_c},any(size(V)~=[n_a,n_z]))
    fprintf('%s: StationaryDist sums to one [T1], this should be zero: %2.8e \n',cmdname{c_c},abs(sum(StationaryDist(:))-1))
    fprintf('%s: V is finite everywhere [T0], this should be zero: %i \n',cmdname{c_c},any(~isfinite(V(:))))
    [~,~,~,mcstatdist]=MarkovChainMoments(z_grid,pi_z);
    zmarg=sum(StationaryDist,1);
    fprintf('%s: z marginal of the agent dist equals the chain stationary dist [T2], this should be small: %2.8e \n',cmdname{c_c},max(abs(zmarg(:)-mcstatdist(:))))
    fprintf('%s: AllStats mean of z equals the chain mean [T2], this should be small: %2.8e \n',cmdname{c_c},abs(AllStats.zvalue.Mean-sum(mcstatdist(:).*z_grid(:))))
    output.cmd(c_c).name=cmdname{c_c};
end

end
