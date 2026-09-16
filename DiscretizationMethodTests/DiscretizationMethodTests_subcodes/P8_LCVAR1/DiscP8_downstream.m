function output=DiscP8_downstream(calib)
% P8: feed this block's outputs into an actual finite-horizon toolkit solve
%
% The first downstream test in the bank with MORE THAN ONE exogenous variable, so it is the first
% that exercises the stacked z_grid_J against n_z as a vector. That pairing is the whole point: the
% command returns sum(znum)-by-J stacked, the solver is told n_z=[znum,znum], and if either side
% disagreed about which block of rows belongs to which variable the age profiles would come back
% attached to the wrong variable - with both variables having similar magnitudes, that is a failure
% that looks like a mild accuracy problem rather than an error.
%
% BOTH THE PLAIN CASE AND THE V_Jplus1 CASE, because the pi_z_J third dimension is J-1. The docs say
% to APPEND a slice when using vfoptions.V_Jplus1 rather than overwrite the last one, and that
% instruction only became correct when the commands moved to J-1 slices (B14). Nothing tested it.

fprintf('\n========== P8: downstream (do the outputs go into a finite-horizon solve) ========== \n')

output=struct();
J=calib.J; M=calib.M; znum=5; % 5^2=25 joint states, enough to be a real solve and small enough to be quick
mewzT=calib.vary.mewz; sigmazT=calib.vary.sigmaz;
fprintf('the truth has variable 1''s variance running %2.4f at age 1 to %2.4f at age %i, and variable 2''s \n',sigmazT(1,1)^2,sigmazT(1,J)^2,J)
fprintf('%2.4f to %2.4f, so the two profiles have different shapes and cannot be confused for each other \n',sigmazT(2,1)^2,sigmazT(2,J)^2)

[z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleVAR1_Tauchen(calib.vary.Mew_J,calib.vary.Rho_J,calib.vary.SigmaSq_J,znum,J,struct());

n_d=0; d_grid=[];
n_a=1; a_grid=1;          % nothing to choose, so the agent's z path IS the chain
n_z=[znum,znum];
Params=struct(); Params.beta=0.96;
DiscountFactorParamNames={'beta'};
ReturnFn=@(aprime,a,z1,z2) z1+z2;
FnsToEvaluate.z1=@(aprime,a,z1,z2) z1;
FnsToEvaluate.z2=@(aprime,a,z1,z2) z2;

vfoptions=struct(); simoptions=struct();
[V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,J,d_grid,a_grid,z_grid_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
jequaloneDist=zeros([n_a,n_z]);
jequaloneDist(1,:)=gather(jequaloneDistz(:))';
Params.ageweights=ones(1,J)/J;
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,{'ageweights'},Policy,n_d,n_a,n_z,J,pi_z_J,Params,simoptions);
AgeStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,J,d_grid,a_grid,z_grid_J,simoptions);

fprintf('\nsolve runs, V has the right number of elements [T0], this should be zero: %i \n',numel(V)~=n_a*prod(n_z)*J)
fprintf('V is finite everywhere [T0], this should be zero: %i \n',any(~isfinite(gather(V(:)))))
fprintf('StationaryDist sums to one [T1], this should be zero: %2.8e \n',abs(sum(gather(StationaryDist(:)))-1))

% The two variables, checked SEPARATELY against their own profiles. If the stacked grid were being
% unpacked the other way round, variable 1's statistics would follow variable 2's truth.
for m_c=1:M
    fn=['z',num2str(m_c)];
    agemean=gather(AgeStats.(fn).Mean); agevar=gather(AgeStats.(fn).Variance);
    [dm,jm]=max(abs(agemean(:)'-mewzT(m_c,:)));
    [dv,jv]=max(abs(agevar(:)'-sigmazT(m_c,:).^2));
    fprintf('variable %i: age-conditional mean, worst age j=%i, error [T2] %2.3e \n',m_c,jm,dm)
    fprintf('variable %i: age-conditional variance, worst age j=%i, error [T2] %2.3e \n',m_c,jv,dv)
    % ...and against the OTHER variable's truth, which must be worse. This is the check that catches
    % a swapped unpacking, and it is not implied by the two lines above.
    other=3-m_c;
    dv_other=max(abs(agevar(:)'-sigmazT(other,:).^2));
    fprintf('variable %i: its variance profile matches ITS OWN truth better than variable %i''s [T2], this should be one: %i \n',m_c,other,dv<dv_other)
    % The age-shift test P6 established: a one-age shift would misreport the whole profile. av is
    % forced to a ROW first - agevar comes back 1-by-J, and mixing a column with a row here makes a
    % (J-1)-by-(J-1) matrix by implicit expansion rather than a scalar, which is what the 2026-09-16
    % run printed 24 times out of this one fprintf.
    % THE AGE-SHIFT TEST RUNS ON THE MEAN, NOT THE VARIANCE, and it took three attempts to get here.
    % The instrument has to be a statistic whose unshifted error is small next to the age-to-age
    % change in the truth, or the test cannot tell a shifted profile from an inaccurate one. At the
    % znum=5 this block runs, the variance carries a genuine truncation bias of 2.3e-02 while the
    % truth changes by about the same amount per age, so the variance fails that requirement - which
    % is why variable 2's variance version kept failing with nothing wrong. The mean since B33 is
    % exact to 8e-17, against age-to-age changes of order 1e-02, so it discriminates by fourteen
    % orders of magnitude. Both maxima are over ages 1..J-1 so they cover the same index set.
    % The mean is only usable where the true mean actually MOVES with age - on a driftless
    % calibration a shifted constant is the same constant and the test would fail with nothing
    % wrong, which is the trap DiscP6_downstream hits with FGP. This block's calibration has drift
    % in both variables by construction, but the instrument is chosen rather than assumed, and which
    % one was used is printed.
    am=agemean(:)'; av=agevar(:)';
    if max(abs(diff(mewzT(m_c,:))))>1e-10
        prof=am; ptruth=mewzT(m_c,:); pname='mean';
    else
        prof=av; ptruth=sigmazT(m_c,:).^2; pname='variance';
    end
    asis_r=max(abs(prof(1:end-1)-ptruth(1:end-1)));
    shifted_r=max(abs(prof(1:end-1)-ptruth(2:end)));
    fprintf('variable %i: the %s profile matches the truth better unshifted than shifted by one age, both over ages 1 to J-1 [T2], this should be one: %i \n',m_c,pname,asis_r<shifted_r)
    fprintf('variable %i: unshifted %2.3e against shifted %2.3e \n',m_c,asis_r,shifted_r)
    output.var(m_c).mean=agemean; output.var(m_c).var=agevar;
end

%% The V_Jplus1 case
% pi_z_J has J-1 slices. docs/ExogenousShocks.md says to APPEND a slice when supplying V_Jplus1, not
% to overwrite the last one - the overwrite instruction was correct only back when the commands
% emitted J slices and the last one was padding. So the appended chain has J slices, and the solve
% has to accept that shape as well.
V_Jplus1=zeros([n_a,n_z]);
pi_z_J_app=cat(3,pi_z_J,pi_z_J(:,:,end)); % append: the age-J transition, reusing the last available one
vfoptions2=struct(); vfoptions2.V_Jplus1=V_Jplus1;
fprintf('\nthe appended transition has J slices rather than J-1 [T0], this should be one: %i \n',size(pi_z_J_app,3)==J)
[V2,~]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,J,d_grid,a_grid,z_grid_J,pi_z_J_app,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('the solve accepts the appended J-slice chain with V_Jplus1 [T0], this should be zero: %i \n',numel(V2)~=n_a*prod(n_z)*J)
fprintf('V is finite everywhere [T0], this should be zero: %i \n',any(~isfinite(gather(V2(:)))))
% With V_Jplus1 all zeros the age-J continuation is worth nothing, which is exactly the terminal
% condition the plain solve imposes, so the two must agree.
fprintf('V_Jplus1=0 gives the same answer as the plain solve [T1], this should be zero: %2.8e \n',max(abs(gather(V2(:))-gather(V(:)))))

output.V=gather(V);

end
