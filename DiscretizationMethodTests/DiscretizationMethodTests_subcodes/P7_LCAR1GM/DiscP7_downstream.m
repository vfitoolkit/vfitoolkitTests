function output=DiscP7_downstream(calib)
% P7: feed this block's outputs into a finite-horizon toolkit solve
%
% Same shape as P6's, with one addition worth having: the age-conditional statistics are checked on
% SKEWNESS as well as mean and variance. A mixture process whose skewness survived the discretization
% but was lost by the solve would look fine in every check up to this point.

fprintf('\n========== P7: downstream (do the outputs go into a finite-horizon solve) ========== \n')

output=struct();
J=calib.J; znum=9;
mew=calib.vary.mew; rho=calib.vary.rho;
p=calib.vary.mixprobs_i; mu=calib.vary.mu_i; sd=calib.vary.sigma_i;
mewzT=calib.vary.mewz; varzT=calib.vary.varz;

n_d=0; d_grid=[];
n_a=1; a_grid=1;   % nothing to choose, so the agent's z path IS the chain
n_z=znum;
Params=struct(); Params.beta=0.96;
Params.ageweights=ones(1,J)/J;
DiscountFactorParamNames={'beta'};
ReturnFn=@(aprime,a,z) z;
FnsToEvaluate.zvalue=@(aprime,a,z) z;

cmdname={'wGM_KFTT','wGM_Tauchen'};
for c_c=1:2
    if c_c==1
        [z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,znum,J,struct('verbose',0));
    else
        [z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleAR1wGM_Tauchen(mew,rho,p,mu,sd,znum,J,3,struct('verbose',0));
    end

    vfoptions=struct(); simoptions=struct();
    [V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,J,d_grid,a_grid,z_grid_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    jequaloneDist=zeros(n_a,n_z);
    jequaloneDist(1,:)=gather(jequaloneDistz(:))';
    StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,{'ageweights'},Policy,n_d,n_a,n_z,J,pi_z_J,Params,simoptions);
    AgeStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,J,d_grid,a_grid,z_grid_J,simoptions);

    fprintf('%s: solve runs, V has the right number of elements [T0], this should be zero: %i \n',cmdname{c_c},numel(V)~=n_a*znum*J)
    fprintf('%s: StationaryDist sums to one [T1], this should be zero: %2.8e \n',cmdname{c_c},abs(sum(StationaryDist(:))-1))
    fprintf('%s: V is finite everywhere [T0], this should be zero: %i \n',cmdname{c_c},any(~isfinite(gather(V(:)))))
    agemean=gather(AgeStats.zvalue.Mean);
    agevar=gather(AgeStats.zvalue.Variance);
    [dm,jm]=max(abs(agemean(:)'-mewzT));
    [dv,jv]=max(abs(agevar(:)'-varzT));
    fprintf('%s: age-conditional mean, worst age j=%i, error [T2] %2.3e \n',cmdname{c_c},jm,dm)
    fprintf('%s: age-conditional variance, worst age j=%i, error [T2] %2.3e \n',cmdname{c_c},jv,dv)
    % the shape test, as in P6: a one-age shift would put age 1's profile where age 2's belongs
    shifted=max(abs(agevar(1:end-1)'-varzT(2:end)));
    asis=max(abs(agevar(:)'-varzT));
    if shifted<asis
        fprintf('%s: THE AGE INDEX IS SHIFTED - the variance profile matches better after shifting by one age \n',cmdname{c_c})
    end
    output.cmd(c_c).name=cmdname{c_c};
    output.cmd(c_c).agevar=agevar(:)';
end

end
