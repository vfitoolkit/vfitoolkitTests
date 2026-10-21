function output=QHDInfHorz_d_z_GI(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% InfHorz quasi-hyperbolic, with d, with z, WITH the grid interpolation layer.
%
% Raws reached:
%   Naive         -> NO quasi-hyperbolic raw at all, as of 2026-10-21. The dispatcher now does the
%                    exponential solve by the standard postGI path (ValueFnIter_InfHorz_Refine_postGI_raw,
%                    via ValueFnIter_InfHorz_GridInterpLayer) and then takes the one quasi-hyperbolic
%                    step with ValueFnIter_Case1_FHorz at N_j=1 and vfoptions.V_Jplus1=Valt, which
%                    gets divide-and-conquer plus the grid interpolation layer with no search
%                    window. ValueFnIter_InfHorz_QuasiHyperbolicN_Refine_postGI_raw was left with
%                    no call site by that change and was DELETED from the toolkit on 2026-10-21.
%   Sophisticated -> ValueFnIter_InfHorz_QuasiHyperbolicS_Refine_postGI_raw
%
% With n_d>0 and gridinterplayer=1 the dispatcher always goes to the refine+postGI raws;
% vfoptions.solnmethod is not consulted at all on this path, so there is no plain variant here.
%
% Policy has four rows: (d, L1, L2, L2flag). The d channel is recovered by looking dstar up at
% the chosen FINE-grid index, which is a separate lookup from the one used for Policyalt. The
% refinement check below verifies it directly: the recorded d must be the argmax over d of the
% return function at the recorded aprime.
%
% These two raws are the only quasi-hyperbolic InfHorz raws that branch on vfoptions.lowmemory,
% so this is the one subcode where the lowmemory check is not a no-op.

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();

N_d=prod(n_d);
N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);

%% Quasi-Hyperbolic Discounting, with grid interpolation
vfoptions.exoticpreferences='QuasiHyperbolic';
vfoptions.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=5;

n2short=vfoptions.ngridinterp;
N_afine=N_a+(N_a-1)*n2short;

beta=Params.beta;
beta0beta=Params.beta0*Params.beta;
avals=repmat(a_grid,1,N_z);
zvals=repmat(z_grid',N_a,1);

%% The exponential-discounting reference solve, same model, same grid interpolation layer
vfoptionsexp=struct(); % deliberately no exoticpreferences field
vfoptionsexp.gridinterplayer=1;
vfoptionsexp.ngridinterp=vfoptions.ngridinterp;
vfoptionsexp.howardssparse=0;
[Vexp,Policyexp]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsexp);

%% ------------------------------------------------------------------ Naive
vfoptions.quasi_hyperbolic='Naive';

vfoptions1=vfoptions;
vfoptions1.howardssparse=0;
simoptions1=simoptions;
simoptions1.gridinterplayer=vfoptions1.gridinterplayer;
simoptions1.ngridinterp=vfoptions1.ngridinterp;
[V1,Policy1,V1alt,Policy1alt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

PolicyVals1=PolicyInd2Val_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,vfoptions1);
PolicyVals1alt=PolicyInd2Val_InfHorz(Policy1alt,n_d,n_a,n_z,d_grid,a_grid,vfoptions1);
fprintf('Naive (refine+GI): PolicyInd2Val ran on both Policy and Policyalt \n')

%% Naive (refine+GI): shapes and index ranges
fprintf('Naive (refine+GI): size(Policy) should be [4 %d %d], is [%s] \n',N_a,N_z,num2str(size(Policy1)))
dpol=gather(reshape(Policy1(1,:,:),[N_a,N_z]));
L1=gather(reshape(Policy1(2,:,:),[N_a,N_z]));
L2=gather(reshape(Policy1(3,:,:),[N_a,N_z]));
L2flag=gather(reshape(Policy1(4,:,:),[N_a,N_z]));
fprintf('Naive (refine+GI): d index in range, this should be one: %d \n',all(dpol(:)>=1) && all(dpol(:)<=N_d))
fprintf('Naive (refine+GI): L1 in [1,N_a-1], this should be one: %d \n',all(L1(:)>=1) && all(L1(:)<=N_a-1))
fprintf('Naive (refine+GI): L2 in [1,ngridinterp+2], this should be one: %d \n',all(L2(:)>=1) && all(L2(:)<=n2short+2))
fprintf('Naive (refine+GI): L2flag in {1,2,3}, this should be one: %d \n',all(ismember(L2flag(:),[1,2,3])))
fprintf('Naive (refine+GI): combined fine index in [1,%d], this should be one: %d \n',N_afine,all((n2short+1)*(L1(:)-1)+L2(:)>=1) && all((n2short+1)*(L1(:)-1)+L2(:)<=N_afine))
dpolalt=gather(reshape(Policy1alt(1,:,:),[N_a,N_z]));
L1alt=gather(reshape(Policy1alt(2,:,:),[N_a,N_z]));
L2alt=gather(reshape(Policy1alt(3,:,:),[N_a,N_z]));
fprintf('Naive (refine+GI): d index in range (alt), this should be one: %d \n',all(dpolalt(:)>=1) && all(dpolalt(:)<=N_d))
fprintf('Naive (refine+GI): L1alt in [1,N_a-1], this should be one: %d \n',all(L1alt(:)>=1) && all(L1alt(:)<=N_a-1))
fprintf('Naive (refine+GI): L2alt in [1,ngridinterp+2], this should be one: %d \n',all(L2alt(:)>=1) && all(L2alt(:)<=n2short+2))

%% Naive (refine+GI): Valt and Policyalt must equal the exponential solve of the same model
fprintf('Naive (refine+GI): Valt vs exponential ValueFnIter with GI, this should be zero: %2.8e \n',max(abs(V1alt(:)-Vexp(:))))
fprintf('Naive (refine+GI): Policyalt vs exponential ValueFnIter with GI, this should be zero: %.3e \n',max(abs(Policy1alt(:)-Policyexp(:))))

%% Naive (refine+GI): the one-step Bellman identity, recomputed here from scratch
%   aprime  = a_grid(L1) + w*(a_grid(L1+1)-a_grid(L1)),  w=(L2-1)/(ngridinterp+1)
%   E[Valt] = (1-w)*EValt(L1,z) + w*EValt(L1+1,z),  EValt on the COARSE grid
%   Vtilde  = u(d,aprime,a,z) + beta0*beta*E[Valt]
% Exact by construction. This pins the discount factor AND the packing of all four Policy rows.
wgt=(L2-1)/(n2short+1);
aprimevals=a_grid(L1)+wgt.*(a_grid(L1+1)-a_grid(L1));
dvals=d_grid(dpol);
c=(1+Params.r)*avals+Params.w*dvals.*zvals-aprimevals;
Fpol=-Inf(N_a,N_z);
okpol=(c>0)&(dvals<1);
Fpol(okpol)=(c(okpol).^(1-Params.sigma)-1)/(1-Params.sigma)+Params.varphi*((1-dvals(okpol)).^(1-Params.eta)-1)/(1-Params.eta);
EValtcoarse=gather(V1alt)*gather(pi_z)';  % N_a-by-N_z: column z is E[Valt(aprime,z')|z] on the coarse aprime grid
EValt=zeros(N_a,N_z);
for z_c=1:N_z
    EValt(:,z_c)=(1-wgt(:,z_c)).*EValtcoarse(L1(:,z_c),z_c)+wgt(:,z_c).*EValtcoarse(L1(:,z_c)+1,z_c);
end
resid=gather(V1)-(Fpol+beta0beta*EValt);
fprintf('Naive (refine+GI): one-step Bellman identity (Vtilde = u + beta0*beta*E[Valt]), this should be zero: %2.8e \n',max(abs(resid(:))))
fprintf('Naive (refine+GI): chosen (d,aprime) always affordable, this should be one: %d \n',all(isfinite(Fpol(:))))
PolicyVals1d=gather(reshape(PolicyVals1(1,:,:),[N_a,N_z]));
PolicyVals1aprime=gather(reshape(PolicyVals1(2,:,:),[N_a,N_z]));
fprintf('Naive (refine+GI): reconstructed d matches PolicyInd2Val, this should be zero: %2.8e \n',max(abs(dvals(:)-PolicyVals1d(:))))
fprintf('Naive (refine+GI): reconstructed aprime matches PolicyInd2Val, this should be zero: %2.8e \n',max(abs(aprimevals(:)-PolicyVals1aprime(:))))
residwrongdiscount=gather(V1)-(Fpol+beta*EValt);
fprintf('Naive (refine+GI): the SAME identity but discounting by beta alone, this should be WELL ABOVE zero: %2.8f \n',max(abs(residwrongdiscount(:))))

%% Naive (refine+GI): the recorded d must be argmax over d at the recorded aprime
% This is the refinement condition, checked directly. It is the dstar lookup at the chosen fine
% index that the check targets: get that lookup wrong and the d recorded is some other state's.
worstdgap=-Inf;
for d_c=1:N_d
    cd=(1+Params.r)*avals+Params.w*d_grid(d_c)*zvals-aprimevals;
    Fd=-Inf(N_a,N_z);
    okd=(cd>0)&(d_grid(d_c)<1);
    Fd(okd)=(cd(okd).^(1-Params.sigma)-1)/(1-Params.sigma)+Params.varphi*((1-d_grid(d_c))^(1-Params.eta)-1)/(1-Params.eta);
    worstdgap=max(worstdgap,max(Fd(:)-Fpol(:)));
end
fprintf('Naive (refine+GI): recorded d is argmax over d at the recorded aprime, this should be zero or below: %2.8e \n',worstdgap)

%% Naive (refine+GI): lowmemory
% Both quasi-hyperbolic refine+postGI raws genuinely branch on lowmemory: lowmemory=1 loops over
% z when building the rough and fine return matrices. Same answer, different construction.
vfoptions1_lm=vfoptions1;
vfoptions1_lm.lowmemory=1;
[V1lm,Policy1lm,V1altlm,Policy1altlm]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_lm);
fprintf('Naive (refine+GI): lowmemory=1, V, this should be zero: %2.8e \n',max(abs(V1(:)-V1lm(:))))
fprintf('Naive (refine+GI): lowmemory=1, Policy, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1lm(:))))
fprintf('Naive (refine+GI): lowmemory=1, Valt, this should be zero: %2.8e \n',max(abs(V1alt(:)-V1altlm(:))))
fprintf('Naive (refine+GI): lowmemory=1, Policyalt, this should be zero: %.3e \n',max(abs(Policy1alt(:)-Policy1altlm(:))))

%% Naive (refine+GI): the maxaprimediff window
% The default with a decision variable is the ordinary exponential discounting one, ceil(n_a/5)
% at this grid size. It bounds the exponential leg only -- the quasi-hyperbolic step is
% unwindowed -- so widening it must not change anything at all, neither leg.
% maxaprimediff=125 with n_a=251 is the WHOLE grid: 2*125+1=251=N_a, and aprimeshifter collapses
% to the single value 126, so the window spans [1,N_a] at every state. That makes this an exact
% unconstrained-search reference. It is also the most memory-hungry solve in the bank: the fine
% return matrix is built at N_d-by-N_aprime-by-N_a-by-N_z before the max over d. If it runs out
% of GPU memory, drop it to ~90 (still far above the default) rather than deleting the check.
vfoptions1_wide=vfoptions1;
vfoptions1_wide.maxaprimediff=125;
[V1w,Policy1w,V1altw,Policy1altw]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_wide);
% maxaprimediff no longer reaches the Naive quasi-hyperbolic step: that step is now a one-period
% FHorz divide-and-conquer + grid-interpolation solve with no search window. So the first two lines
% are no longer a window-WIDTH check on the quasi-hyperbolic step; they are the stronger statement
% that the step is window-INDEPENDENT, and an exact identity rather than something we hope is
% small. If either moves, maxaprimediff has leaked back into the quasi-hyperbolic step. Without
% this relabelling the block would be a vacuous pass. The third line is still a genuine window
% check, for the exponential leg, which does still window.
fprintf('Naive (refine+GI): quasi-hyperbolic step is window-independent, Vtilde, this should be zero: %2.8e \n',max(abs(V1(:)-V1w(:))))
fprintf('Naive (refine+GI): quasi-hyperbolic step is window-independent, Policy, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1w(:))))
fprintf('Naive (refine+GI): exponential leg under a wider window, Valt, this should be zero: %2.8e \n',max(abs(V1alt(:)-V1altw(:))))
% The width measurements below are now DIAGNOSTICS ONLY on this path: nothing Naive depends on
% maxaprimediff any more. They are kept because they are the only measurement anywhere of the
% quasi-hyperbolic-vs-exponential policy gap, which is the evidence the Sophisticated no-default
% refusal message cites for recommending ceil(n_a/(beta0*4)).
%
% NO LONGER CENSORED: aprimeshifter's clamp applied to the windowed quasi-hyperbolic step, and that
% step is now unwindowed, so this is the gap itself rather than a lower bound on it. It should come
% out EQUAL to the whole-grid line below. That is not an independent check -- identical Policy was
% already established above -- it is reported so the two lines cannot be misread as disagreeing.
fprintf('Naive (refine+GI): max |L1 - L1alt| at the DEFAULT window, DIAGNOSTIC (no longer censored, the step is unwindowed): %d \n',max(abs(L1(:)-L1alt(:))))
% The whole-grid run, which never had a clamp on either leg. This is the width a window WOULD need.
L1w=gather(reshape(Policy1w(2,:,:),[N_a,N_z]));
L1altw=gather(reshape(Policy1altw(2,:,:),[N_a,N_z]));
fprintf('Naive (refine+GI): max |L1 - L1alt| on the WHOLE-GRID run, UNCENSORED, this is the width the window needs: %d \n',max(abs(L1w(:)-L1altw(:))))
% SECOND beta0 -- see the matching block in QHDInfHorz_nod_z_GI.m for why one calibration cannot
% decide whether the default rule should key off n_a or off beta0.
Params_b2=Params;
Params_b2.beta0=0.80;
[~,Policy1w2,~,Policy1altw2]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_b2,DiscountFactorParamNames,[],vfoptions1_wide);
L1w2=gather(reshape(Policy1w2(2,:,:),[N_a,N_z]));
L1altw2=gather(reshape(Policy1altw2(2,:,:),[N_a,N_z]));
fprintf('Naive (refine+GI): required width at beta0=0.80 (against beta0=%2.2f above), UNCENSORED: %d \n',Params.beta0,max(abs(L1w2(:)-L1altw2(:))))
clear Policy1w2 Policy1altw2

%% Naive (refine+GI): howards
vfoptions1_noH=vfoptions1;
vfoptions1_noH.howards=0;
[V1noH,Policy1noH]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_noH);
fprintf('Naive (refine+GI): howards=0 (pure VFI), V, this should be zero: %2.8e \n',max(abs(V1(:)-V1noH(:))))
fprintf('Naive (refine+GI): howards=0 (pure VFI), Policy, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1noH(:))))

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
fprintf('Naive (refine+GI): ValueFnFromPolicy V (Vtilde), this should be zero to convergence tolerance: %2.8e \n',max(abs(VfromPolicy(:)-V1(:))))
fprintf('Naive (refine+GI): ValueFnFromPolicy Valt (V_std), this should be zero to convergence tolerance: %2.8e \n',max(abs(ValtfromPolicy(:)-V1alt(:))))

%% Naive (refine+GI): anti-vacuity
fprintf('Naive (refine+GI): quasi-hyperbolic matters, max|Vtilde-V_std| should be well above zero: %2.8f \n',max(abs(V1(:)-V1alt(:))))
fprintf('Naive (refine+GI): quasi-hyperbolic matters, max|Policy-Policyalt| should be well above zero: %2.8f \n',max(abs(Policy1(:)-Policy1alt(:))))
fprintf('Naive (refine+GI): the d channel is not degenerate, number of distinct d chosen should be above one: %d \n',numel(unique(dpol(:))))
fprintf('Naive (refine+GI): the L2 channel is not degenerate, number of distinct L2 chosen should be above one: %d \n',numel(unique(L2(:))))

clear V1lm V1w V1noH Policy1lm Policy1w Policy1noH V1altlm V1altw Policy1altlm

%% ---------------------------------------------------------- Sophisticated
% Reminder: maxiter is forced to 1000 and there is no Howards in any Sophisticated raw.
vfoptions2=vfoptions;
vfoptions2.quasi_hyperbolic='Sophisticated';
% Sophisticated + grid interpolation has NO default maxaprimediff (the dispatcher errors unless
% it is set). Use the value the error message recommends, so the bank exercises it.
vfoptions2.maxaprimediff=ceil(N_a/(Params.beta0*4));
vfoptions2.howardssparse=0;
simoptions2=simoptions;
simoptions2.gridinterplayer=vfoptions2.gridinterplayer;
simoptions2.ngridinterp=vfoptions2.ngridinterp;
[V2,Policy2,V2alt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
% V2=Vhat, V2alt=Vunderbar

PolicyVals2=PolicyInd2Val_InfHorz(Policy2,n_d,n_a,n_z,d_grid,a_grid,vfoptions2);
fprintf('Sophisticated (refine+GI): PolicyInd2Val ran \n')

fprintf('Sophisticated (refine+GI): size(Policy) should be [4 %d %d], is [%s] \n',N_a,N_z,num2str(size(Policy2)))
dpols=gather(reshape(Policy2(1,:,:),[N_a,N_z]));
L1s=gather(reshape(Policy2(2,:,:),[N_a,N_z]));
L2s=gather(reshape(Policy2(3,:,:),[N_a,N_z]));
L2flags=gather(reshape(Policy2(4,:,:),[N_a,N_z]));
fprintf('Sophisticated (refine+GI): d index in range, this should be one: %d \n',all(dpols(:)>=1) && all(dpols(:)<=N_d))
fprintf('Sophisticated (refine+GI): L1 in [1,N_a-1], this should be one: %d \n',all(L1s(:)>=1) && all(L1s(:)<=N_a-1))
fprintf('Sophisticated (refine+GI): L2 in [1,ngridinterp+2], this should be one: %d \n',all(L2s(:)>=1) && all(L2s(:)<=n2short+2))
fprintf('Sophisticated (refine+GI): L2flag in {1,2,3}, this should be one: %d \n',all(ismember(L2flags(:),[1,2,3])))
fprintf('Sophisticated (refine+GI): Vhat is finite everywhere, this should be one: %d \n',gather(all(isfinite(V2(:)))))
fprintf('Sophisticated (refine+GI): Vunderbar is finite everywhere, this should be one: %d \n',gather(all(isfinite(V2alt(:)))))

%% Sophisticated (refine+GI): the two value functions at the single equilibrium policy
wgts=(L2s-1)/(n2short+1);
aprimevals2=a_grid(L1s)+wgts.*(a_grid(L1s+1)-a_grid(L1s));
dvals2=d_grid(dpols);
c2=(1+Params.r)*avals+Params.w*dvals2.*zvals-aprimevals2;
Fpol2=-Inf(N_a,N_z);
okpol2=(c2>0)&(dvals2<1);
Fpol2(okpol2)=(c2(okpol2).^(1-Params.sigma)-1)/(1-Params.sigma)+Params.varphi*((1-dvals2(okpol2)).^(1-Params.eta)-1)/(1-Params.eta);
EVunderbarcoarse=gather(V2alt)*gather(pi_z)';
EVunderbar=zeros(N_a,N_z);
for z_c=1:N_z
    EVunderbar(:,z_c)=(1-wgts(:,z_c)).*EVunderbarcoarse(L1s(:,z_c),z_c)+wgts(:,z_c).*EVunderbarcoarse(L1s(:,z_c)+1,z_c);
end
residhat=gather(V2)-(Fpol2+beta0beta*EVunderbar);
residunderbar=gather(V2alt)-(Fpol2+beta*EVunderbar);
% The CONVERGENCE MEASUREMENT, not checks -- see the matching block in QHDInfHorz_nod_z.m for why
% these must not carry check wording.
fprintf('Sophisticated (refine+GI): Vhat residual against u + beta0*beta*E[Vunderbar], DIAGNOSTIC not a check: %2.8e \n',max(abs(residhat(:))))
fprintf('Sophisticated (refine+GI): Vunderbar residual against u + beta*E[Vunderbar], THE CONVERGENCE DIAGNOSTIC, not a check: %2.8e \n',max(abs(residunderbar(:))))
fprintf('Sophisticated (refine+GI): chosen (d,aprime) always affordable, this should be one: %d \n',all(isfinite(Fpol2(:))))
% The identity above says Vhat-Vunderbar = (beta0*beta-beta)*E[Vunderbar], and beta0*beta<beta,
% so Vhat exceeds Vunderbar exactly where E[Vunderbar] is NEGATIVE -- not everywhere. With
% sigma=2 the flow utility is u=1-1/c, which is POSITIVE for c>1, so V is not signed in this
% calibration and 'Vhat>Vunderbar everywhere' is false in the nod models and true in the d
% models (where the leisure term drags V negative). The share below is the correct check.
signok=(gather(V2)-gather(V2alt)>0)==(EVunderbar<0);
fprintf('Sophisticated (refine+GI): Vhat exceeds Vunderbar exactly where E[Vunderbar]<0, share of states should be one: %2.4f \n',mean(signok(:)))

%% Sophisticated (refine+GI): lowmemory
vfoptions2_lm=vfoptions2;
vfoptions2_lm.lowmemory=1;
[V2lm,Policy2lm,V2altlm]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2_lm);
fprintf('Sophisticated (refine+GI): lowmemory=1, Vhat, this should be zero: %2.8e \n',max(abs(V2(:)-V2lm(:))))
fprintf('Sophisticated (refine+GI): lowmemory=1, Vunderbar, this should be zero: %2.8e \n',max(abs(V2alt(:)-V2altlm(:))))
fprintf('Sophisticated (refine+GI): lowmemory=1, Policy, this should be zero: %.3e \n',max(abs(Policy2(:)-Policy2lm(:))))

%% Sophisticated (refine+GI): ValueFnFromPolicy_InfHorz
% Sophisticated needs no Policyalt: Policy alone determines both Vhat and Vunderbar.
[VfromPolicyS,ValtfromPolicyS]=ValueFnFromPolicy_InfHorz(Policy2,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions2);
fprintf('Sophisticated (refine+GI): ValueFnFromPolicy V (Vhat) vs the solver, CONVERGENCE-DEPENDENT (this measures solver non-convergence, it is not a fault): %2.8e \n',max(abs(VfromPolicyS(:)-V2(:))))
fprintf('Sophisticated (refine+GI): ValueFnFromPolicy Valt (Vunderbar) vs the solver, CONVERGENCE-DEPENDENT: %2.8e \n',max(abs(ValtfromPolicyS(:)-V2alt(:))))

%% Sophisticated (refine+GI): anti-vacuity
fprintf('Sophisticated (refine+GI): quasi-hyperbolic matters, max|Vhat-Vunderbar| should be well above zero: %2.8f \n',max(abs(V2(:)-V2alt(:))))
fprintf('Sophisticated (refine+GI): quasi-hyperbolic matters, max|Policy-Policyexp| should be well above zero: %2.8f \n',max(abs(Policy2(:)-Policyexp(:))))

%% Naive vs Sophisticated
fprintf('Naive vs Sophisticated (refine+GI): max|V_Naive - V_Sophisticated| should be well above zero: %2.8f \n',max(abs(V1(:)-V2(:))))

%% Figures
zmid=ceil(N_z/2);

fig=figure(figure_c);
subplot(1,3,1); plot(a_grid,squeeze(PolicyVals1(2,:,zmid)), a_grid,squeeze(PolicyVals1alt(2,:,zmid)), a_grid,a_grid,'k:')
title('Naive, refine+GI: aprime policy (mid z)'); xlabel('a'); legend('quasi-hyperbolic','exponential','45 degree','Location','northwest')
subplot(1,3,2); plot(a_grid,squeeze(PolicyVals1(1,:,zmid)), a_grid,squeeze(PolicyVals1alt(1,:,zmid)))
title('Naive, refine+GI: d policy (mid z)'); xlabel('a'); legend('quasi-hyperbolic','exponential')
subplot(1,3,3); plot(a_grid,V1(:,zmid), a_grid,V1alt(:,zmid))
title('Naive, refine+GI: value fn (mid z)'); xlabel('a'); legend('Vtilde','V\_std','Location','southeast')

fig=figure(100+figure_c);
subplot(1,3,1); plot(a_grid,squeeze(PolicyVals2(2,:,zmid)), a_grid,squeeze(PolicyVals1alt(2,:,zmid)), a_grid,a_grid,'k:')
title('Sophisticated, refine+GI: aprime policy (mid z)'); xlabel('a'); legend('quasi-hyperbolic','exponential','45 degree','Location','northwest')
subplot(1,3,2); plot(a_grid,squeeze(PolicyVals2(1,:,zmid)), a_grid,squeeze(PolicyVals1alt(1,:,zmid)))
title('Sophisticated, refine+GI: d policy (mid z)'); xlabel('a'); legend('quasi-hyperbolic','exponential')
subplot(1,3,3); plot(a_grid,V2(:,zmid), a_grid,V2alt(:,zmid))
title('Sophisticated, refine+GI: value fn (mid z)'); xlabel('a'); legend('Vhat','Vunderbar','Location','southeast')

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
