function output=DiscP10_downstream(calib)
% P10: feed the AR(p) outputs into an actual toolkit solve
%
% The point of this block is not that the solver runs - P5 already established that a stacked
% multi-dimensional z_grid with n_z a vector goes through. It is that the LAG STRUCTURE survives the
% trip. Dimension 2 of the state is dimension 1 one period earlier, so in the stationary agent
% distribution the two must have the same marginal, and the model's own statistics must say so.
% A command that wired the shift up wrongly would produce a perfectly good-looking solve whose
% second state variable was not the lag of the first, and nothing else in the bank would notice.
%
% The return function deliberately uses BOTH z1 and z2, so a mis-wired lag changes the policy rather
% than just sitting in the exogenous block unread.

fprintf('\n========== P10: downstream (do the outputs go into a solve) ========== \n')

output=struct();
Rho=calib.AR2.Rho; p=length(Rho);
mew=calib.mew; sigma=calib.sigma;
T=calib.AR2.gauss;
znum=7;

n_d=0; d_grid=[];
n_a=51; a_grid=linspace(0,20,n_a)';
n_z=znum*ones(1,p);
Params=struct(); Params.beta=0.94; Params.r=0.04;
DiscountFactorParamNames={'beta'};
ReturnFn=@(aprime,a,z1,z2,r) log(max((1+r)*a+exp(z1)+0.25*exp(z2)-aprime,10^(-10)));
FnsToEvaluate.z1=@(aprime,a,z1,z2) z1;
FnsToEvaluate.z2=@(aprime,a,z1,z2) z2;
FnsToEvaluate.assets=@(aprime,a,z1,z2) a;

cmdname={'ARp_FarmerToda','ARpwGM_FarmerToda'};
for c_c=1:2
    fo=struct(); fo.verbose=0;
    if c_c==1
        [z_grid,pi_z]=discretizeARp_FarmerToda(mew,Rho,sigma,znum,fo);
        Tc=T;
    else
        [z_grid,pi_z]=discretizeARpwGM_FarmerToda(mew,Rho,calib.mixprobs_i,calib.mu_i,calib.sigma_i,znum,fo);
        Tc=calib.AR2.gm;
    end

    vfoptions=struct(); simoptions=struct();
    [V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params,[]);
    AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

    fprintf('\n%s: solve runs, V has the right number of elements [T0], this should be zero: %i \n',cmdname{c_c},numel(V)~=n_a*prod(n_z))
    fprintf('%s: V is finite everywhere [T0], this should be zero: %i \n',cmdname{c_c},any(~isfinite(gather(V(:)))))
    fprintf('%s: StationaryDist sums to one [T1], this should be zero: %2.8e \n',cmdname{c_c},abs(sum(gather(StationaryDist(:)))-1))

    m1=gather(AllStats.z1.Mean); v1=gather(AllStats.z1.Variance);
    m2=gather(AllStats.z2.Mean); v2=gather(AllStats.z2.Variance);
    % THE LAG IDENTITY. z2 is z1 one period earlier, so in the stationary distribution the two are
    % the same random variable and every marginal moment must agree. This is the check the block
    % exists for, and it is an identity rather than an accuracy claim - it holds however coarse the
    % grid is, because it is a property of the chain and not of how well the chain approximates.
    fprintf('%s: z2 has the same stationary mean as z1, since it IS z1 lagged [T1], this should be zero: %2.8e \n',cmdname{c_c},abs(m2-m1))
    fprintf('%s: and the same variance [T1], this should be zero: %2.8e \n',cmdname{c_c},abs(v2-v1))
    % ...and both against the analytic truth, which is a T2 accuracy statement and a separate thing.
    fprintf('%s: E(z) %+2.6f against truth %+2.6f, error [T2] %2.3e \n',cmdname{c_c},m1,Tc.Ez,abs(m1-Tc.Ez))
    fprintf('%s: Var(z) %2.6f against truth %2.6f, error [T2] %2.3e \n',cmdname{c_c},v1,Tc.var,abs(v1-Tc.var))
    % DID THE MODEL ACTUALLY USE THE EXOGENOUS STATES? The first version of this asked whether the
    % asset distribution was non-degenerate, which turned out to test the calibration rather than
    % the toolkit: with beta*(1+r)=0.977 and an income standard deviation of about 0.13 the agent
    % is impatient and sits on the borrowing constraint, so mean assets are zero and the variance is
    % numerical noise - it passed for one command and failed for the other on the run of 2026-09-17,
    % decided by which side of 1e-10 the noise fell. The question worth asking is whether the two
    % exogenous states reach the value function at all, and that is answered directly: V must vary
    % along BOTH z dimensions. It holds whatever the asset policy does.
    Vr=reshape(gather(V),[n_a,znum,znum]);
    varyz1=max(abs(Vr-Vr(:,1,:)),[],'all');
    varyz2=max(abs(Vr-Vr(:,:,1)),[],'all');
    fprintf('%s: V varies along z1, so the first exogenous state reaches the return function [T0], this should be one: %i \n',cmdname{c_c},varyz1>1e-08)
    fprintf('%s: V varies along z2, so the LAG does too - it is not along for the ride [T0], this should be one: %i \n',cmdname{c_c},varyz2>1e-08)
    fprintf('%s: mean assets %2.4f (the agent is impatient here, so this sits at the constraint) \n',cmdname{c_c},gather(AllStats.assets.Mean))

    output.cmd(c_c).name=cmdname{c_c};
    output.cmd(c_c).m1=m1; output.cmd(c_c).v1=v1; output.cmd(c_c).m2=m2; output.cmd(c_c).v2=v2;
end

end
