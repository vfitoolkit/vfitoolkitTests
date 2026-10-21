function output=QHDInfHorz_nod_z_GI(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% InfHorz quasi-hyperbolic, without d, with z, WITH the grid interpolation layer.
%
% Raws reached:
%   Naive         -> NO quasi-hyperbolic raw at all, as of 2026-10-21. The dispatcher now does the
%                    exponential solve by the standard postGI path and then takes the one
%                    quasi-hyperbolic step with ValueFnIter_Case1_FHorz at N_j=1 and
%                    vfoptions.V_Jplus1=Valt, which gets divide-and-conquer plus the grid
%                    interpolation layer with no search window.
%                    ValueFnIter_InfHorz_QuasiHyperbolicN_postGI_nod_raw was left with no call
%                    site by that change and was DELETED from the toolkit on 2026-10-21.
%   Sophisticated -> ValueFnIter_InfHorz_QuasiHyperbolicS_postGI_nod_raw
%
% The quasi-hyperbolic grid interpolation layer is always postGI: the dispatcher sets the postGI
% defaults itself and never reads vfoptions.preGI, so preGI=1 is silently ignored here.
%
% Policy has three rows: (L1, L2, L2flag). L1 is the lower coarse-grid index, L2 the index into
% the ngridinterp+2 points spanning [L1,L1+1], and L2flag the -Inf-avoidance flag. Both the
% aprime value and the interpolated continuation are LINEAR in the weight (L2-1)/(ngridinterp+1),
% which is what makes the exact Bellman identity below reproducible by hand.
%
% The maxaprimediff sensitivity check used to be the one that mattered most on this path, because
% the QH step took its fine-grid argmax inside a window of +-vfoptions.maxaprimediff coarse points
% around the EXPONENTIAL policy, and a quasi-hyperbolic policy that drifted further than that was
% silently pinned at the window edge. Since 2026-10-21 the Naive QH step has no window at all, so
% that check has become the stronger claim that the step is window-INDEPENDENT. The window still
% bounds the exponential leg, where vfoptions.maxaprimediff now takes the ordinary exponential
% discounting default (mirrored from ValueFnIter_InfHorz_GridInterpLayer).

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
% Do the current setup
n_d=0; d_grid=[];

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);

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

%% The exponential-discounting reference solve, same model, same grid interpolation layer.
% howardssparse is forced to 0 on BOTH so that the exponential leg of the quasi-hyperbolic raw
% and the reference use the same Howards implementation (the quasi-hyperbolic postGI raws have
% no sparse variant, while the standard GI default is howardssparse=1).
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
fprintf('Naive (GI): PolicyInd2Val ran on both Policy and Policyalt \n')

%% Naive (GI): shapes, index ranges, and the L2/L2flag channels
fprintf('Naive (GI): size(Policy) should be [3 %d %d], is [%s] \n',N_a,N_z,num2str(size(Policy1)))
L1=gather(reshape(Policy1(1,:,:),[N_a,N_z]));
L2=gather(reshape(Policy1(2,:,:),[N_a,N_z]));
L2flag=gather(reshape(Policy1(3,:,:),[N_a,N_z]));
fprintf('Naive (GI): L1 in [1,N_a-1], this should be one: %d \n',all(L1(:)>=1) && all(L1(:)<=N_a-1))
fprintf('Naive (GI): L2 in [1,ngridinterp+2], this should be one: %d \n',all(L2(:)>=1) && all(L2(:)<=n2short+2))
fprintf('Naive (GI): L2flag in {1,2,3}, this should be one: %d \n',all(ismember(L2flag(:),[1,2,3])))
fprintf('Naive (GI): combined fine index in [1,%d], this should be one: %d \n',N_afine,all((n2short+1)*(L1(:)-1)+L2(:)>=1) && all((n2short+1)*(L1(:)-1)+L2(:)<=N_afine))
L1alt=gather(reshape(Policy1alt(1,:,:),[N_a,N_z]));
L2alt=gather(reshape(Policy1alt(2,:,:),[N_a,N_z]));
L2flagalt=gather(reshape(Policy1alt(3,:,:),[N_a,N_z]));
fprintf('Naive (GI): L1alt in [1,N_a-1], this should be one: %d \n',all(L1alt(:)>=1) && all(L1alt(:)<=N_a-1))
fprintf('Naive (GI): L2alt in [1,ngridinterp+2], this should be one: %d \n',all(L2alt(:)>=1) && all(L2alt(:)<=n2short+2))
fprintf('Naive (GI): L2flagalt in {1,2,3}, this should be one: %d \n',all(ismember(L2flagalt(:),[1,2,3])))

%% Naive (GI): Valt and Policyalt must equal the exponential solve of the same model
fprintf('Naive (GI): Valt vs exponential ValueFnIter with GI, this should be zero: %2.8e \n',max(abs(V1alt(:)-Vexp(:))))
fprintf('Naive (GI): Policyalt vs exponential ValueFnIter with GI, this should be zero: %.3e \n',max(abs(Policy1alt(:)-Policyexp(:))))

%% Naive (GI): the one-step Bellman identity, recomputed here from scratch
% Both the chosen aprime and the interpolated continuation are linear in the within-segment
% weight w=(L2-1)/(ngridinterp+1), so the whole quasi-hyperbolic step can be rebuilt by hand:
%   aprime  = a_grid(L1) + w*(a_grid(L1+1)-a_grid(L1))
%   E[Valt] = (1-w)*EValt(L1,z) + w*EValt(L1+1,z),  EValt on the COARSE grid
%   Vtilde  = u(aprime,a,z) + beta0*beta*E[Valt]
% Exact by construction: the raw takes exactly one fine-grid argmax step off the returned Valt.
wgt=(L2-1)/(n2short+1);
aprimevals=a_grid(L1)+wgt.*(a_grid(L1+1)-a_grid(L1));
c=(1+Params.r)*avals+Params.w*zvals-aprimevals;
Fpol=-Inf(N_a,N_z);
Fpol(c>0)=(c(c>0).^(1-Params.sigma)-1)/(1-Params.sigma);
EValtcoarse=gather(V1alt)*gather(pi_z)';  % N_a-by-N_z: column z is E[Valt(aprime,z')|z] over the coarse aprime grid
EValt=zeros(N_a,N_z);
for z_c=1:N_z
    EValt(:,z_c)=(1-wgt(:,z_c)).*EValtcoarse(L1(:,z_c),z_c)+wgt(:,z_c).*EValtcoarse(L1(:,z_c)+1,z_c);
end
resid=gather(V1)-(Fpol+beta0beta*EValt);
fprintf('Naive (GI): one-step Bellman identity (Vtilde = u + beta0*beta*E[Valt]), this should be zero: %2.8e \n',max(abs(resid(:))))
fprintf('Naive (GI): chosen aprime always affordable, this should be one: %d \n',all(isfinite(Fpol(:))))
PolicyVals1aprime=gather(reshape(PolicyVals1(1,:,:),[N_a,N_z]));
fprintf('Naive (GI): reconstructed aprime matches PolicyInd2Val, this should be zero: %2.8e \n',max(abs(aprimevals(:)-PolicyVals1aprime(:))))
residwrongdiscount=gather(V1)-(Fpol+beta*EValt);
fprintf('Naive (GI): the SAME identity but discounting by beta alone, this should be WELL ABOVE zero: %2.8f \n',max(abs(residwrongdiscount(:))))

%% Naive (GI): the maxaprimediff window
% The default here is the ordinary exponential discounting default for a sole-a1prime choice,
% which is a small constant either side of the previous iterate's own argmax. It bounds the
% exponential leg only -- the quasi-hyperbolic step is unwindowed -- so widening it must not
% change anything at all, neither leg.
% maxaprimediff=125 with n_a=251 is the WHOLE grid: 2*125+1=251=N_a, and aprimeshifter collapses
% to the single value 126, so the window spans [1,N_a] at every state. That makes this an exact
% unconstrained-search reference rather than merely a wider window.
vfoptions1_wide=vfoptions1;
vfoptions1_wide.maxaprimediff=125;
[V1w,Policy1w,V1altw,Policy1altw]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_wide);
% maxaprimediff no longer reaches the Naive quasi-hyperbolic step at all: that step is now a
% one-period FHorz divide-and-conquer + grid-interpolation solve with no search window. So the
% first two lines are no longer a window-WIDTH check on the quasi-hyperbolic step; they are the
% stronger statement that the step is window-INDEPENDENT, and an exact identity rather than
% something we hope comes out small. If either ever moves, maxaprimediff has leaked back into the
% quasi-hyperbolic step. Without this relabelling the block would be a vacuous pass: it would read
% as zero while testing nothing about the step it names.
% The third line is still a genuine window check, for the EXPONENTIAL leg, which does still window.
fprintf('Naive (GI): quasi-hyperbolic step is window-independent, Vtilde, this should be zero: %2.8e \n',max(abs(V1(:)-V1w(:))))
fprintf('Naive (GI): quasi-hyperbolic step is window-independent, Policy, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1w(:))))
fprintf('Naive (GI): exponential leg under a wider window, Valt, this should be zero: %2.8e \n',max(abs(V1alt(:)-V1altw(:))))
% How far the quasi-hyperbolic policy actually sits from the exponential one, in coarse points.
% This is no longer a constraint to check against the default: the quasi-hyperbolic step is
% unwindowed, so the gap may exceed the default window freely without pinning anything.
% The three width measurements below are now DIAGNOSTICS ONLY on this path: nothing Naive depends
% on maxaprimediff any more. They are kept because they are the only measurement anywhere of the
% quasi-hyperbolic-vs-exponential policy gap, and that is the evidence the Sophisticated
% no-default refusal message cites for recommending ceil(n_a/(beta0*4)). Delete them and that
% formula becomes unsourced.
%
% NO LONGER CENSORED. In the windowed era aprimeshifter clamped the quasi-hyperbolic step, so this
% statistic could not exceed the default and had to be read as a lower bound. The step is now
% unwindowed, so nothing clamps L1 and this is the gap itself. It should therefore come out EQUAL
% to the whole-grid line below -- which is not an independent check, because the two runs were
% already shown to give identical Policy above; it is reported so the two lines cannot be misread
% as disagreeing.
fprintf('Naive (GI): max |L1 - L1alt| at the DEFAULT window, DIAGNOSTIC (no longer censored, the step is unwindowed): %d \n',max(abs(L1(:)-L1alt(:))))
% The whole-grid run, which never had a clamp on either leg. This is the width a window WOULD have
% to span, and so the number the Sophisticated refusal message's formula is sourced from.
L1w=gather(reshape(Policy1w(1,:,:),[N_a,N_z]));
L1altw=gather(reshape(Policy1altw(1,:,:),[N_a,N_z]));
fprintf('Naive (GI): max |L1 - L1alt| on the WHOLE-GRID run, UNCENSORED, this is the width the window needs: %d \n',max(abs(L1w(:)-L1altw(:))))
% SECOND beta0, same uncensored measurement. The default rule keys off n_a alone, but the width the
% window has to span is the quasi-hyperbolic-vs-exponential POLICY GAP, and that widens as beta0
% falls. One calibration cannot tell the two apart. If the two numbers below are close, an
% n_a-proportional rule is right and only its constant is in question; if the lower beta0 needs a
% lot more, the rule has to key off beta0 too and no single n_a fraction will do.
Params_b2=Params;
Params_b2.beta0=0.80;
[~,Policy1w2,~,Policy1altw2]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_b2,DiscountFactorParamNames,[],vfoptions1_wide);
L1w2=gather(reshape(Policy1w2(1,:,:),[N_a,N_z]));
L1altw2=gather(reshape(Policy1altw2(1,:,:),[N_a,N_z]));
fprintf('Naive (GI): required width at beta0=0.80 (against beta0=%2.2f above), UNCENSORED: %d \n',Params.beta0,max(abs(L1w2(:)-L1altw2(:))))
clear Policy1w2 Policy1altw2

%% Naive (GI): ngridinterp
% More interpolation points must not move the coarse part of the answer much, and must move the
% value function towards the fine-grid solve.
vfoptions1_ng=vfoptions1;
vfoptions1_ng.ngridinterp=11;
[V1ng,Policy1ng]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_ng);
fprintf('Naive (GI): ngridinterp=11 vs 5, Vtilde, this should be small (both approximate the same continuous problem): %2.8e \n',max(abs(V1(:)-V1ng(:))))

%% Naive (GI): does the grid interpolation layer actually buy accuracy?
% Solve the same model by plain pure discretization on the big (1001 point) grid, which contains
% the coarse grid at every 10th point, and compare. GI on the coarse grid should sit closer to
% the big-grid answer than plain discretization on the coarse grid does.
vfoptions1_plain=vfoptions;
vfoptions1_plain.gridinterplayer=0;
vfoptions1_plain=rmfield(vfoptions1_plain,'ngridinterp');
[V1c,Policy1c]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_plain);
[V1big,Policy1big]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_plain);
V1bigcoarse=V1big(1:((n_a_big-1)/(n_a-1)):end,:);
fprintf('Naive: |V(coarse grid, no GI) - V(big grid)| : %2.8f \n',max(abs(gather(V1c(:))-gather(V1bigcoarse(:)))))
fprintf('Naive: |V(coarse grid, with GI) - V(big grid)|, should be the SMALLER of the two: %2.8f \n',max(abs(gather(V1(:))-gather(V1bigcoarse(:)))))

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
fprintf('Naive (GI): ValueFnFromPolicy V (Vtilde), this should be zero to convergence tolerance: %2.8e \n',max(abs(VfromPolicy(:)-V1(:))))
fprintf('Naive (GI): ValueFnFromPolicy Valt (V_std), this should be zero to convergence tolerance: %2.8e \n',max(abs(ValtfromPolicy(:)-V1alt(:))))

%% Naive (GI): anti-vacuity
fprintf('Naive (GI): quasi-hyperbolic matters, max|Vtilde-V_std| should be well above zero: %2.8f \n',max(abs(V1(:)-V1alt(:))))
fprintf('Naive (GI): quasi-hyperbolic matters, max|aprime - aprimealt| should be well above zero: %2.8f \n',max(abs(PolicyVals1(:)-PolicyVals1alt(:))))
fprintf('Naive (GI): the L2 channel is not degenerate, number of distinct L2 chosen should be above one: %d \n',numel(unique(L2(:))))

clear V1w V1ng V1c V1big Policy1w Policy1ng Policy1c Policy1big V1altw Policy1altw

%% ---------------------------------------------------------- Sophisticated
% Reminder: maxiter is forced to 1000 and there is no Howards in any Sophisticated raw.
vfoptions2=vfoptions;
vfoptions2.quasi_hyperbolic='Sophisticated';
% Sophisticated + grid interpolation has NO default maxaprimediff: the dispatcher errors unless it
% is set, because the width depends on beta0 and a too-narrow window fails silently. The bank sets
% it to the value the error message recommends, so that this run exercises the recommendation --
% if ceil(n_a/(beta0*4)) is not enough, the checks below are what will say so.
vfoptions2.maxaprimediff=ceil(N_a/(Params.beta0*4));
vfoptions2.howardssparse=0;
simoptions2=simoptions;
simoptions2.gridinterplayer=vfoptions2.gridinterplayer;
simoptions2.ngridinterp=vfoptions2.ngridinterp;
[V2,Policy2,V2alt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
% V2=Vhat, V2alt=Vunderbar

PolicyVals2=PolicyInd2Val_InfHorz(Policy2,n_d,n_a,n_z,d_grid,a_grid,vfoptions2);
fprintf('Sophisticated (GI): PolicyInd2Val ran \n')

fprintf('Sophisticated (GI): size(Policy) should be [3 %d %d], is [%s] \n',N_a,N_z,num2str(size(Policy2)))
L1s=gather(reshape(Policy2(1,:,:),[N_a,N_z]));
L2s=gather(reshape(Policy2(2,:,:),[N_a,N_z]));
L2flags=gather(reshape(Policy2(3,:,:),[N_a,N_z]));
fprintf('Sophisticated (GI): L1 in [1,N_a-1], this should be one: %d \n',all(L1s(:)>=1) && all(L1s(:)<=N_a-1))
fprintf('Sophisticated (GI): L2 in [1,ngridinterp+2], this should be one: %d \n',all(L2s(:)>=1) && all(L2s(:)<=n2short+2))
fprintf('Sophisticated (GI): L2flag in {1,2,3}, this should be one: %d \n',all(ismember(L2flags(:),[1,2,3])))
fprintf('Sophisticated (GI): Vhat is finite everywhere, this should be one: %d \n',gather(all(isfinite(V2(:)))))
fprintf('Sophisticated (GI): Vunderbar is finite everywhere, this should be one: %d \n',gather(all(isfinite(V2alt(:)))))

%% Sophisticated (GI): the two value functions at the single equilibrium policy
wgts=(L2s-1)/(n2short+1);
aprimevals2=a_grid(L1s)+wgts.*(a_grid(L1s+1)-a_grid(L1s));
c2=(1+Params.r)*avals+Params.w*zvals-aprimevals2;
Fpol2=-Inf(N_a,N_z);
Fpol2(c2>0)=(c2(c2>0).^(1-Params.sigma)-1)/(1-Params.sigma);
EVunderbarcoarse=gather(V2alt)*gather(pi_z)';
EVunderbar=zeros(N_a,N_z);
for z_c=1:N_z
    EVunderbar(:,z_c)=(1-wgts(:,z_c)).*EVunderbarcoarse(L1s(:,z_c),z_c)+wgts(:,z_c).*EVunderbarcoarse(L1s(:,z_c)+1,z_c);
end
residhat=gather(V2)-(Fpol2+beta0beta*EVunderbar);
residunderbar=gather(V2alt)-(Fpol2+beta*EVunderbar);
% The CONVERGENCE MEASUREMENT, not checks -- see the matching block in QHDInfHorz_nod_z.m for why
% these must not carry check wording.
fprintf('Sophisticated (GI): Vhat residual against u + beta0*beta*E[Vunderbar], DIAGNOSTIC not a check: %2.8e \n',max(abs(residhat(:))))
fprintf('Sophisticated (GI): Vunderbar residual against u + beta*E[Vunderbar], THE CONVERGENCE DIAGNOSTIC, not a check: %2.8e \n',max(abs(residunderbar(:))))
fprintf('Sophisticated (GI): chosen aprime always affordable, this should be one: %d \n',all(isfinite(Fpol2(:))))
% The identity above says Vhat-Vunderbar = (beta0*beta-beta)*E[Vunderbar], and beta0*beta<beta,
% so Vhat exceeds Vunderbar exactly where E[Vunderbar] is NEGATIVE -- not everywhere. With
% sigma=2 the flow utility is u=1-1/c, which is POSITIVE for c>1, so V is not signed in this
% calibration and 'Vhat>Vunderbar everywhere' is false in the nod models and true in the d
% models (where the leisure term drags V negative). The share below is the correct check.
signok=(gather(V2)-gather(V2alt)>0)==(EVunderbar<0);
fprintf('Sophisticated (GI): Vhat exceeds Vunderbar exactly where E[Vunderbar]<0, share of states should be one: %2.4f \n',mean(signok(:)))

%% Sophisticated (GI): maxaprimediff
vfoptions2_wide=vfoptions2;
vfoptions2_wide.maxaprimediff=125;
[V2w,Policy2w,V2altw]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2_wide);
% NOT a window-width check for Sophisticated, and not claimed as one. At beta0<1 the continuation
% the quasi-hyperbolic objective interpolates is Vunderbar -- one criterion (beta) evaluated at
% another criterion's (beta0*beta) argmax -- so it is kinked, and LINEAR interpolation of a kinked
% function collapses the fine argmax onto an interval endpoint however wide the window is. Widening
% maxaprimediff therefore moves the answer here for a reason that has nothing to do with the window,
% and asserting zero would be asserting something the solution concept does not deliver. The Naive
% lines in the sibling subcodes ARE genuine window-width checks; these are the measurement of a
% known, documented limitation. At beta0=1 Vunderbar IS Vhat and nothing collapses.
fprintf('Sophisticated (GI): maxaprimediff vs the whole grid, Vhat, DIAGNOSTIC (kinked-Vunderbar collapse, not window width): %2.8e \n',max(abs(V2(:)-V2w(:))))
fprintf('Sophisticated (GI): maxaprimediff vs the whole grid, Vunderbar, DIAGNOSTIC: %2.8e \n',max(abs(V2alt(:)-V2altw(:))))
fprintf('Sophisticated (GI): maxaprimediff vs the whole grid, Policy, DIAGNOSTIC: %.3e \n',max(abs(Policy2(:)-Policy2w(:))))

%% Sophisticated (GI): ValueFnFromPolicy_InfHorz
% Sophisticated needs no Policyalt: Policy alone determines both Vhat and Vunderbar.
[VfromPolicyS,ValtfromPolicyS]=ValueFnFromPolicy_InfHorz(Policy2,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions2);
fprintf('Sophisticated (GI): ValueFnFromPolicy V (Vhat) vs the solver, CONVERGENCE-DEPENDENT (this measures solver non-convergence, it is not a fault): %2.8e \n',max(abs(VfromPolicyS(:)-V2(:))))
fprintf('Sophisticated (GI): ValueFnFromPolicy Valt (Vunderbar) vs the solver, CONVERGENCE-DEPENDENT: %2.8e \n',max(abs(ValtfromPolicyS(:)-V2alt(:))))

%% Sophisticated (GI): anti-vacuity
fprintf('Sophisticated (GI): quasi-hyperbolic matters, max|Vhat-Vunderbar| should be well above zero: %2.8f \n',max(abs(V2(:)-V2alt(:))))
fprintf('Sophisticated (GI): quasi-hyperbolic matters, max|Policy-Policyexp| should be well above zero: %2.8f \n',max(abs(Policy2(:)-Policyexp(:))))

%% Naive vs Sophisticated
fprintf('Naive vs Sophisticated (GI): max|V_Naive - V_Sophisticated| should be well above zero: %2.8f \n',max(abs(V1(:)-V2(:))))

%% Figures
zmid=ceil(N_z/2);

fig=figure(figure_c);
subplot(1,2,1); plot(a_grid,squeeze(PolicyVals1(1,:,zmid)), a_grid,squeeze(PolicyVals1alt(1,:,zmid)), a_grid,a_grid,'k:')
title('Naive, GI: aprime policy (mid z)'); xlabel('a'); legend('quasi-hyperbolic','exponential','45 degree','Location','northwest')
subplot(1,2,2); plot(a_grid,V1(:,zmid), a_grid,V1alt(:,zmid))
title('Naive, GI: value fn (mid z)'); xlabel('a'); legend('Vtilde','V\_std','Location','southeast')

fig=figure(100+figure_c);
subplot(1,2,1); plot(a_grid,squeeze(PolicyVals2(1,:,zmid)), a_grid,squeeze(PolicyVals1alt(1,:,zmid)), a_grid,a_grid,'k:')
title('Sophisticated, GI: aprime policy (mid z)'); xlabel('a'); legend('quasi-hyperbolic','exponential','45 degree','Location','northwest')
subplot(1,2,2); plot(a_grid,V2(:,zmid), a_grid,V2alt(:,zmid))
title('Sophisticated, GI: value fn (mid z)'); xlabel('a'); legend('Vhat','Vunderbar','Location','southeast')

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
