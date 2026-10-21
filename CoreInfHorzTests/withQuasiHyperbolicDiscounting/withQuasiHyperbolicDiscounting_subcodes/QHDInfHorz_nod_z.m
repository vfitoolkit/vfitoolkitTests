function output=QHDInfHorz_nod_z(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% InfHorz quasi-hyperbolic, without d, with z, pure discretization (no grid interpolation layer).
%
% Raws reached:
%   Naive         -> ValueFnIter_InfHorz_QuasiHyperbolicN_nod_raw   (after ValueFnIter_InfHorz_nod_raw for V_std)
%   Sophisticated -> ValueFnIter_InfHorz_QuasiHyperbolicS_nod_raw
%
% The centrepiece here is the one-step Bellman identity, recomputed in this file from Policy,
% Valt, a_grid, z_grid and pi_z without calling the toolkit at all. For Naive it is exact by
% construction (the raw takes a single argmax step off the returned Valt), so it pins the
% discount factor that multiplies the continuation. That is the check that catches a Naive
% raw discounting by beta0 rather than beta0*beta.
%
% Note: ValueFnFromPolicy_InfHorz is not yet available as a quasi-hyperbolic oracle. It errors on
% any vfoptions carrying an 'exoticpreferences' field, so these hand-written identities are what
% does that job for now. The ValueFnFromPolicy calls are written test-first further down, and
% become a second, independent oracle once the toolkit support is built.

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
% Do the current setup
n_d=0; d_grid=[];

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);

%% Quasi-Hyperbolic Discounting
vfoptions.exoticpreferences='QuasiHyperbolic';
vfoptions.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;

% Useful for the identity checks below
beta=Params.beta;
beta0beta=Params.beta0*Params.beta;
avals=repmat(a_grid,1,N_z);
zvals=repmat(z_grid',N_a,1);

%% The exponential-discounting reference solve for this same model.
% For Naive, Valt/Policyalt ARE the exponential solution (the QH step is taken off it), so this
% must match exactly. For Sophisticated it is just a reference for the figure and for the
% anti-vacuity check.
vfoptionsexp=struct(); % deliberately no exoticpreferences field
[Vexp,Policyexp]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsexp);

%% ------------------------------------------------------------------ Naive
vfoptions.quasi_hyperbolic='Naive';

vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1,V1alt,Policy1alt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

PolicyVals1=PolicyInd2Val_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,vfoptions1);
PolicyVals1alt=PolicyInd2Val_InfHorz(Policy1alt,n_d,n_a,n_z,d_grid,a_grid,vfoptions1);
fprintf('Naive: PolicyInd2Val ran on both Policy and Policyalt \n')

% Shapes and index ranges
fprintf('Naive: size(V) should be [%d %d], is [%s] \n',N_a,N_z,num2str(size(V1)))
fprintf('Naive: size(Policy) should be [1 %d %d], is [%s] \n',N_a,N_z,num2str(size(Policy1)))
fprintf('Naive: aprime index in range, this should be one: %d \n',gather(all(Policy1(:)>=1) && all(Policy1(:)<=N_a)))
fprintf('Naive: aprimealt index in range, this should be one: %d \n',gather(all(Policy1alt(:)>=1) && all(Policy1alt(:)<=N_a)))
fprintf('Naive: V is finite everywhere, this should be one: %d \n',gather(all(isfinite(V1(:)))))
fprintf('Naive: Valt is finite everywhere, this should be one: %d \n',gather(all(isfinite(V1alt(:)))))

%% Naive: Valt and Policyalt must equal the plain exponential-discounting solution
% The Naive raw is handed the converged exponential V_std and takes one step off it, so these
% are the very same numbers produced by the exponential solve of the same model.
fprintf('Naive: Valt vs plain exponential ValueFnIter, this should be zero: %.3e \n',max(abs(V1alt(:)-Vexp(:))))
fprintf('Naive: Policyalt vs plain exponential ValueFnIter, this should be zero: %.3e \n',max(abs(Policy1alt(:)-Policyexp(:))))

%% Naive: the one-step Bellman identity, recomputed here from scratch
% Vtilde(a,z) = u(aprime(a,z),a,z) + beta0*beta*E[Valt(aprime(a,z),z')|z]
% Exact by construction: the raw takes exactly one argmax step off the returned Valt.
aprimeindex=gather(reshape(Policy1,[N_a,N_z]));
aprimevals=a_grid(aprimeindex);
c=(1+Params.r)*avals+Params.w*zvals-aprimevals;
Fpol=-Inf(N_a,N_z);
Fpol(c>0)=(c(c>0).^(1-Params.sigma)-1)/(1-Params.sigma);
EValt=zeros(N_a,N_z,'gpuArray');
for z_c=1:N_z
    EValt(:,z_c)=V1alt(aprimeindex(:,z_c),:)*pi_z(z_c,:)';
end
resid=V1-(Fpol+beta0beta*EValt);
fprintf('Naive: one-step Bellman identity (Vtilde = u + beta0*beta*E[Valt]), this should be zero: %2.8e \n',max(abs(resid(:))))
fprintf('Naive: chosen aprime always affordable, this should be one: %d \n',gather(all(isfinite(Fpol(:)))))

% Same identity, but with beta in place of beta0*beta, so that the diary records how far the
% two discount factors are apart. If a raw ever discounts by the wrong one, THIS line is the
% one that goes to zero and the line above is the one that does not.
residwrongdiscount=V1-(Fpol+beta*EValt);
fprintf('Naive: the SAME identity but discounting by beta alone, this should be WELL ABOVE zero: %2.8f \n',max(abs(residwrongdiscount(:))))
residwrongdiscount2=V1-(Fpol+Params.beta0*EValt);
fprintf('Naive: the SAME identity but discounting by beta0 alone, this should be WELL ABOVE zero: %2.8f \n',max(abs(residwrongdiscount2(:))))

%% Naive: Valt satisfies the standard Bellman equation at Policyalt (to convergence tolerance)
aprimealtindex=gather(reshape(Policy1alt,[N_a,N_z]));
aprimealtvals=a_grid(aprimealtindex);
calt=(1+Params.r)*avals+Params.w*zvals-aprimealtvals;
Fpolalt=-Inf(N_a,N_z);
Fpolalt(calt>0)=(calt(calt>0).^(1-Params.sigma)-1)/(1-Params.sigma);
EValtalt=zeros(N_a,N_z,'gpuArray');
for z_c=1:N_z
    EValtalt(:,z_c)=V1alt(aprimealtindex(:,z_c),:)*pi_z(z_c,:)';
end
residalt=V1alt-(Fpolalt+beta*EValtalt);
fprintf('Naive: standard Bellman identity for Valt (V_std = u + beta*E[V_std]), this should be zero to convergence tolerance: %2.8e \n',max(abs(residalt(:))))

%% Naive: optimality of the QH step
% Vtilde is the max, so it must be at least as large as the same object evaluated at Policyalt.
fprintf('Naive: Vtilde >= value of the QH objective at Policyalt, this should be one: %d \n',gather(all(V1(:)>=(Fpolalt(:)+beta0beta*EValtalt(:))-10^(-10))))

%% ValueFnFromPolicy_InfHorz -- a second, independent reconstruction oracle
% Implemented toolkit-side on 2026-08-30: ValueFnFromPolicy_InfHorz is now declared varargout and
% dispatches to ValueFnFromPolicy_InfHorz_QuasiHyperbolic ahead of the grid-interpolation
% hand-off, mirroring ValueFnFromPolicy_FHorz. It returns two outputs:
%   Naive:         [V, Valt] = [Vtilde, exponential-discounter value at Policyalt].
%                  Requires vfoptions.Policyalt, because Naive stores Policy at the beta0*beta
%                  argmax while Valt is the exponential value at the std argmax, so Policy alone
%                  cannot reconstruct both. Policyalt is the 4th output of ValueFnIter for Naive.
%   Sophisticated: [V, Valt] = [Vhat, Vunderbar]. Policy alone suffices; no Policyalt.
% These are now real checks, and they are INDEPENDENT of the hand-written Bellman identities
% above: the identities evaluate the recursion once at the returned V, whereas this re-solves the
% policy-evaluation fixed point from scratch and only then takes the quasi-hyperbolic step.
%
% WHAT TO EXPECT, and it differs by preference type:
%   Naive         -- near zero (~1e-8). The solver's Valt is a converged exponential-discounting
%                    solution, and ValueFnFromPolicy re-solves the same policy evaluation, so the
%                    two agree up to their respective convergence tolerances.
%   Sophisticated -- NOT near zero, and that is not a fault in either command. ValueFnFromPolicy
%                    solves the Vunderbar policy-evaluation fixed point properly, whereas the
%                    solver stopped at the forced maxiter=1000 with a large currdist. The gap is
%                    therefore roughly currdist/(1-beta), and it MEASURES the solver's
%                    non-convergence. It is the sharpest diagnostic of it in the bank: it says how
%                    far the returned Vunderbar is from the true value of its own policy.
vfoptions_vfp=vfoptions1;
vfoptions_vfp.Policyalt=Policy1alt; % Naive quasi-hyperbolic ValueFnFromPolicy needs the exponential-discounter argmax
[VfromPolicy,ValtfromPolicy]=ValueFnFromPolicy_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions_vfp);
fprintf('Naive: ValueFnFromPolicy V (Vtilde), this should be zero to convergence tolerance: %2.8e \n',max(abs(VfromPolicy(:)-V1(:))))
fprintf('Naive: ValueFnFromPolicy Valt (V_std), this should be zero to convergence tolerance: %2.8e \n',max(abs(ValtfromPolicy(:)-V1alt(:))))

%% Naive: anti-vacuity
% If beta0 were calibrated so that the QH and exponential policies coincided everywhere, every
% check above would pass on a degenerate model. These lines are what make that visible.
fprintf('Naive: quasi-hyperbolic matters, max|Policy-Policyalt| should be well above zero: %2.8f \n',max(abs(Policy1(:)-Policy1alt(:))))
fprintf('Naive: quasi-hyperbolic matters, max|Vtilde-V_std| should be well above zero: %2.8f \n',max(abs(V1(:)-V1alt(:))))
fprintf('Naive: policy is not pinned at the borrowing constraint, share of states with aprime>1 should be well above zero: %2.4f \n',gather(mean(Policy1(:)>1)))

%% Naive: Howards improvement is only an accelerator, so turning it off must not change anything
% (it accelerates the V_std leg; the QH step itself has no Howards)
vfoptions1_noH=vfoptions1;
vfoptions1_noH.howards=0;
[V1noH,Policy1noH,V1altnoH,Policy1altnoH]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_noH);
fprintf('Naive: howards=0 (pure VFI), this should be zero: %.3e \n',max(abs(V1(:)-V1noH(:))))
fprintf('Naive: howards=0 (pure VFI), this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1noH(:))))
fprintf('Naive: howards=0 (Valt), this should be zero: %.3e \n',max(abs(V1alt(:)-V1altnoH(:))))
fprintf('Naive: howards=0 (Policyalt), this should be zero: %.3e \n',max(abs(Policy1alt(:)-Policy1altnoH(:))))

%% Naive: howardssparse=1 uses a different Howards implementation for the V_std leg
% The QH plain-discretization branch of the dispatcher calls ValueFnIter_InfHorz_nod_raw
% directly and never looks at vfoptions.howardssparse or vfoptions.howardsgreedy, so this is
% currently a no-op. It is here so that a future dispatcher that DOES honour them is checked.
vfoptions1_sp=vfoptions1;
vfoptions1_sp.howardssparse=1;
[V1sp,Policy1sp]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_sp);
fprintf('Naive: howardssparse=1, this should be zero: %.3e \n',max(abs(V1(:)-V1sp(:))))
fprintf('Naive: howardssparse=1, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1sp(:))))

%% ---------------------------------------------------------- Sophisticated
% Note: the dispatcher forces vfoptions.maxiter=1000 for Sophisticated, and none of the four
% Sophisticated raw here (S_nod_raw) runs Howards on Vunderbar as of 2026-09-01; before that it
% was pure value function iteration capped at 999
% steps. InfHorz Sophisticated is documented as struggling to converge, so nothing below
% asserts tight convergence; the Vunderbar Bellman residual printed here IS the convergence
% diagnostic (the forced vfoptions.verbose=1 prints nothing, as no QH raw reads verbose).
vfoptions2=vfoptions;
vfoptions2.quasi_hyperbolic='Sophisticated';
simoptions2=simoptions;
[V2,Policy2,V2alt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
% V2=Vhat, V2alt=Vunderbar, Policy2=the equilibrium choice. Sophisticated returns no Policyalt.

PolicyVals2=PolicyInd2Val_InfHorz(Policy2,n_d,n_a,n_z,d_grid,a_grid,vfoptions2);
fprintf('Sophisticated: PolicyInd2Val ran \n')

fprintf('Sophisticated: size(V) should be [%d %d], is [%s] \n',N_a,N_z,num2str(size(V2)))
fprintf('Sophisticated: size(Policy) should be [1 %d %d], is [%s] \n',N_a,N_z,num2str(size(Policy2)))
fprintf('Sophisticated: aprime index in range, this should be one: %d \n',gather(all(Policy2(:)>=1) && all(Policy2(:)<=N_a)))
fprintf('Sophisticated: Vhat is finite everywhere, this should be one: %d \n',gather(all(isfinite(V2(:)))))
fprintf('Sophisticated: Vunderbar is finite everywhere, this should be one: %d \n',gather(all(isfinite(V2alt(:)))))

%% Sophisticated: the two value functions at the single equilibrium policy
% Vhat      = u + beta0*beta*E[Vunderbar]
% Vunderbar = u + beta     *E[Vunderbar]
% Both hold only up to the fixed-point residual of the Vunderbar iteration, so these are the
% convergence diagnostic as much as they are a correctness check.
aprimeindex2=gather(reshape(Policy2,[N_a,N_z]));
aprimevals2=a_grid(aprimeindex2);
c2=(1+Params.r)*avals+Params.w*zvals-aprimevals2;
Fpol2=-Inf(N_a,N_z);
Fpol2(c2>0)=(c2(c2>0).^(1-Params.sigma)-1)/(1-Params.sigma);
EVunderbar=zeros(N_a,N_z,'gpuArray');
for z_c=1:N_z
    EVunderbar(:,z_c)=V2alt(aprimeindex2(:,z_c),:)*pi_z(z_c,:)';
end
residhat=V2-(Fpol2+beta0beta*EVunderbar);
residunderbar=V2alt-(Fpol2+beta*EVunderbar);
% These two are the CONVERGENCE MEASUREMENT, not checks, and must not carry check wording. The
% Sophisticated iteration stops at the forced maxiter=1000 with a residual around 1e-2, so a diary
% summary that grades a 'should be near zero' phrase at its own bar will fail them every run for
% something that is documented behaviour. Saying 'this IS the convergence diagnostic' in prose does
% not help: the summary reads the wording, not the parenthetical.
sophresid=gather(max(abs(residunderbar(:))));
fprintf('Sophisticated: Vhat residual against u + beta0*beta*E[Vunderbar], DIAGNOSTIC not a check: %2.8e \n',max(abs(residhat(:))))
fprintf('Sophisticated: Vunderbar residual against u + beta*E[Vunderbar], THE CONVERGENCE DIAGNOSTIC, not a check: %2.8e \n',sophresid)
fprintf('Sophisticated: chosen aprime always affordable, this should be one: %d \n',gather(all(isfinite(Fpol2(:)))))

% The gap between the two is EXACTLY (beta0*beta-beta)*E[Vunderbar] at the fixed point, so unlike
% the two residuals above this one IS an identity -- but only once Vunderbar has converged, because
% the solver forms Vhat against the PREVIOUS Vunderbar. It is therefore claimed as an identity only
% when the convergence diagnostic above says it has earned that, and reported as a diagnostic
% otherwise. Asserting it unconditionally claims something that is not true of a run stopped early.
residgap=(V2-V2alt)-(beta0beta-beta)*EVunderbar;
if sophresid<10^(-6)
    fprintf('Sophisticated: converged, so the Vhat-Vunderbar gap identity is claimed as one \n')
    fprintf('Sophisticated: (Vhat-Vunderbar) minus (beta0*beta-beta)*E[Vunderbar], this should be zero: %2.8e \n',max(abs(residgap(:))))
else
    fprintf('Sophisticated: did NOT converge, so the Vhat-Vunderbar gap identity holds only to that residual \n')
    fprintf('Sophisticated: (Vhat-Vunderbar) minus (beta0*beta-beta)*E[Vunderbar], CONVERGENCE-DEPENDENT: %2.8e \n',max(abs(residgap(:))))
end
% The identity above says Vhat-Vunderbar = (beta0*beta-beta)*E[Vunderbar], and beta0*beta<beta,
% so Vhat exceeds Vunderbar exactly where E[Vunderbar] is NEGATIVE -- not everywhere. With
% sigma=2 the flow utility is u=1-1/c, which is POSITIVE for c>1, so V is not signed in this
% calibration and 'Vhat>Vunderbar everywhere' is false in the nod models and true in the d
% models (where the leisure term drags V negative). The share below is the correct check.
signok=(gather(V2)-gather(V2alt)>0)==(EVunderbar<0);
fprintf('Sophisticated: Vhat exceeds Vunderbar exactly where E[Vunderbar]<0, share of states should be one: %2.4f \n',mean(signok(:)))

%% Sophisticated: ValueFnFromPolicy_InfHorz
% Sophisticated needs no Policyalt: Policy alone determines both Vhat and Vunderbar.
[VfromPolicyS,ValtfromPolicyS]=ValueFnFromPolicy_InfHorz(Policy2,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions2);
fprintf('Sophisticated: ValueFnFromPolicy V (Vhat) vs the solver, CONVERGENCE-DEPENDENT (this measures solver non-convergence, it is not a fault): %2.8e \n',max(abs(VfromPolicyS(:)-V2(:))))
fprintf('Sophisticated: ValueFnFromPolicy Valt (Vunderbar) vs the solver, CONVERGENCE-DEPENDENT: %2.8e \n',max(abs(ValtfromPolicyS(:)-V2alt(:))))

%% Sophisticated: anti-vacuity
fprintf('Sophisticated: quasi-hyperbolic matters, max|Policy-Policyexp| should be well above zero: %2.8f \n',max(abs(Policy2(:)-Policyexp(:))))
fprintf('Sophisticated: quasi-hyperbolic matters, max|Vhat-Vunderbar| should be well above zero: %2.8f \n',max(abs(V2(:)-V2alt(:))))
fprintf('Sophisticated: policy is not pinned at the borrowing constraint, share of states with aprime>1 should be well above zero: %2.4f \n',gather(mean(Policy2(:)>1)))

%% Sophisticated vs Naive: the two solution concepts must differ
fprintf('Naive vs Sophisticated: max|V_Naive - V_Sophisticated| should be well above zero: %2.8f \n',max(abs(V1(:)-V2(:))))

%% Figures
zmid=ceil(N_z/2);

fig=figure(figure_c);
subplot(1,2,1); plot(a_grid,squeeze(PolicyVals1(1,:,zmid)), a_grid,squeeze(PolicyVals1alt(1,:,zmid)), a_grid,a_grid,'k:')
title('Naive: aprime policy (mid z)'); xlabel('a'); legend('quasi-hyperbolic','exponential','45 degree','Location','northwest')
subplot(1,2,2); plot(a_grid,V1(:,zmid), a_grid,V1alt(:,zmid))
title('Naive: value fn (mid z)'); xlabel('a'); legend('Vtilde','V\_std','Location','southeast')

fig=figure(100+figure_c);
subplot(1,2,1); plot(a_grid,squeeze(PolicyVals2(1,:,zmid)), a_grid,squeeze(PolicyVals1alt(1,:,zmid)), a_grid,a_grid,'k:')
title('Sophisticated: aprime policy (mid z)'); xlabel('a'); legend('quasi-hyperbolic','exponential','45 degree','Location','northwest')
subplot(1,2,2); plot(a_grid,V2(:,zmid), a_grid,V2alt(:,zmid))
title('Sophisticated: value fn (mid z)'); xlabel('a'); legend('Vhat','Vunderbar','Location','southeast')

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
