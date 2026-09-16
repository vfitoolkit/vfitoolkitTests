function output=DiscP6_downstream(calib)
% P6: feed this block's outputs into an actual finite-horizon toolkit solve
%
% The first downstream test in the bank that is finite-horizon, and the first that uses
% jequaloneDistz for its actual purpose - the age-1 distribution of the agent. Everything earlier
% only ever checked that it summed to one.
%
% Two things can only fail here. The pi_z_J third dimension is J-1, so a solver that expected J
% slices would either error or silently read the wrong transition at the last age. And the age
% indexing runs one way: pi_z_J(:,:,j) is the transition OUT OF age j. If it were being read as the
% transition INTO age j, the whole age profile would shift by one - which is invisible in a
% stationary block and is exactly what the rising variance profile here can catch.

fprintf('\n========== P6: downstream (do the outputs go into a finite-horizon solve) ========== \n')

output=struct();
J=calib.J; znum=9;
mew=calib.vary.mew; rho=calib.vary.rho; sigma=calib.vary.sigma;
sigmazT=calib.vary.sigmaz; mewzT=calib.vary.mewz;
fprintf('the truth has a variance of z rising from %2.4f at age 1 to %2.4f at age %i - a shifted age \n',sigmazT(1)^2,sigmazT(J)^2,J)
fprintf('index would report the profile of the wrong age, so the shape is the test, not just the level \n')

n_d=0; d_grid=[];
n_a=1; a_grid=1;   % nothing to choose, so the agent's z path IS the chain
n_z=znum;
Params=struct(); Params.beta=0.96;
DiscountFactorParamNames={'beta'};
ReturnFn=@(aprime,a,z) z;
FnsToEvaluate.zvalue=@(aprime,a,z) z;

cmdname={'LCAR1_KFTT','LCAR1_FGP','LCAR1_FGPTauchen'};
for c_c=1:3
    if c_c==1
        [z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleAR1_KFTT(mew,rho,sigma,znum,J,struct());
    elseif c_c==2
        [z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleAR1_FellaGallipoliPan(rho,sigma,znum,J,struct());
    else
        [z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(mew,rho,sigma,znum,J,struct());
    end

    vfoptions=struct(); simoptions=struct();
    [V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,J,d_grid,a_grid,z_grid_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    jequaloneDist=zeros(n_a,n_z);
    jequaloneDist(1,:)=gather(jequaloneDistz(:))';
    Params.ageweights=ones(1,J)/J; % StationaryDist_FHorz_Case1 takes the NAME of the age-weight parameter, not its values
    StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,{'ageweights'},Policy,n_d,n_a,n_z,J,pi_z_J,Params,simoptions);
    AgeStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,J,d_grid,a_grid,z_grid_J,simoptions);

    fprintf('%s: solve runs, V has the right number of elements [T0], this should be zero: %i \n',cmdname{c_c},numel(V)~=n_a*znum*J)
    fprintf('%s: StationaryDist sums to one [T1], this should be zero: %2.8e \n',cmdname{c_c},abs(sum(StationaryDist(:))-1))
    fprintf('%s: V is finite everywhere [T0], this should be zero: %i \n',cmdname{c_c},any(~isfinite(gather(V(:)))))
    agemean=gather(AgeStats.zvalue.Mean);
    agevar=gather(AgeStats.zvalue.Variance);
    if c_c==2
        mtruth=zeros(1,J); % FGP is driftless: it takes no mew
    else
        mtruth=mewzT;
    end
    [dm,jm]=max(abs(agemean(:)'-mtruth));
    [dv,jv]=max(abs(agevar(:)'-sigmazT.^2));
    fprintf('%s: age-conditional mean, worst age j=%i, error [T2] %2.3e \n',cmdname{c_c},jm,dm)
    fprintf('%s: age-conditional variance, worst age j=%i, error [T2] %2.3e \n',cmdname{c_c},jv,dv)
    % The shape test. A one-age shift would put the age-1 variance where age 2's belongs, and the
    % profile is steep at the start, so the two are far apart.
    % av is forced to a ROW before any of this. agevar comes back 1-by-J, so agevar(1:end-1)'
    % was a COLUMN and sigmazT(2:end).^2 a row - implicit expansion made a (J-1)-by-(J-1) matrix,
    % max() of which is a row vector, and "if shifted<asis" on a vector needs every element true.
    % The test could therefore never fire, whatever the profile did. Found by P8, which prints the
    % same comparison unconditionally and got 24 lines of output out of one fprintf.
    % ON THE MEAN, NOT THE VARIANCE, for the reason DiscP8_downstream sets out: the instrument has to
    % be a statistic whose unshifted error is small next to the age-to-age change in the truth. On
    % the 2026-09-16 run these three commands' mean errors were 1.0e-06, 2.2e-16 and 2.8e-16 while
    % their variance errors were orders larger, so the mean is the discriminating one here too. Both
    % maxima are over ages 1..J-1 so they cover the same index set - taking one over all J ages and
    % the other over J-1 lets the shifted version win by dropping whichever age the worst sits at,
    % which is how P8's copy failed with nothing wrong.
    % ...but the mean is only usable where the true mean actually MOVES with age. FGP takes no mew
    % and is driftless, so its true mean is zero at every age and a shifted zero vector is the same
    % zero vector - the comparison would be between two identical numbers and the test would fail
    % with nothing wrong. So the instrument is chosen rather than assumed, and which one was used is
    % printed, because a test that silently switches what it measures is worse than no test.
    am=agemean(:)'; av=agevar(:)';
    if max(abs(diff(mtruth)))>1e-10
        prof=am; ptruth=mtruth; pname='mean';
    else
        prof=av; ptruth=sigmazT.^2; pname='variance';
    end
    asis=max(abs(prof(1:end-1)-ptruth(1:end-1)));
    shifted=max(abs(prof(1:end-1)-ptruth(2:end)));
    fprintf('%s: the %s profile matches the truth better unshifted than shifted by one age, both over ages 1 to J-1 [T2], this should be one: %i \n',cmdname{c_c},pname,asis<shifted)
    fprintf('%s: unshifted %2.3e against shifted %2.3e \n',cmdname{c_c},asis,shifted)
    output.cmd(c_c).name=cmdname{c_c};
    output.cmd(c_c).agevar=agevar(:)';
    output.cmd(c_c).shifted=(shifted<asis);
end

end
