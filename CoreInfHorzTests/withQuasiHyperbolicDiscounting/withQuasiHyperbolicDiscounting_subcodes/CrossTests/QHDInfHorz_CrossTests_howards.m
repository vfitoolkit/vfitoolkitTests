function output=QHDInfHorz_CrossTests_howards(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test: Howards on vs off, for SOPHISTICATED quasi-hyperbolic.
%
% This is the regression test for the Howards fix of 2026-09-01
% (InfHorz_QHSophisticated_Howards_proposal.md steps 2-4). Before that fix, the Howards block in
% ValueFnIter_InfHorz_QuasiHyperbolicS_raw iterated Vunderbar from Vhat rather than from
% Vunderbar, which composes to
%     Vunderbar = F + beta*E[F] + beta*beta0beta*E^2[Vunderbar]
% -- a different fixed point, with an extra flow payoff and a two-period discount. The same block
% sat commented out in ...S_nod_raw and is now live there too.
%
% WHY THIS CANNOT BE A PLAIN "on == off" CHECK
% Howards is an accelerator, so for a CONVERGENT iteration on and off must agree to tolerance and
% the test would be trivial to write. InfHorz Sophisticated does not converge: it stops at the
% forced maxiter=1000 with currdist around 1e-2. Two runs that stop early on different trajectories
% legitimately land in different places, so max|V_on - V_off| is not a zero-check and asserting it
% as one would produce a number that means nothing.
%
% WHAT IS TESTED INSTEAD: SELF-CONSISTENCY, which does not require convergence.
% Whatever the solve returns, the pair (Policy, Vunderbar) should satisfy
%     Vunderbar = the true policy-evaluation value of Policy
% and ValueFnFromPolicy_InfHorz computes exactly that, independently and to its own tolerance. The
% old bug drove Vunderbar toward a DIFFERENT fixed point during the Howards phase, so it biased
% that residual in a way that Howards-off could not. Comparing the residual with Howards on
% against the residual with Howards off therefore isolates the accelerator from the cycling:
%   after the fix -> the two residuals are comparable (both limited by cycling, not by Howards)
%   with the bug  -> the Howards-on residual is materially the worse of the two
% Each residual is computed against that run's OWN policy, so the two runs having different
% policies does not matter.
%
% Paths covered, one per Sophisticated raw that now runs Howards:
%   (a) nod, pure discretization        -> ValueFnIter_InfHorz_QuasiHyperbolicS_nod_raw
%   (b) d,   refinement (the default)   -> ValueFnIter_InfHorz_QuasiHyperbolicS_nod_raw
%   (c) d,   pure discretization        -> ValueFnIter_InfHorz_QuasiHyperbolicS_raw  <- the live bug
%   (d) nod, grid interpolation        -> ValueFnIter_InfHorz_QuasiHyperbolicS_postGI_nod_raw
%   (e) d,   refine + grid interp      -> ValueFnIter_InfHorz_QuasiHyperbolicS_Refine_postGI_raw
% (d) and (e) were added 2026-10-21 when step 5 gave those two raws Howards; before that they had
% none and were explicitly excluded here with the note that followed:
% The two postGI Sophisticated raws are NOT covered: they have no Howards (deferred, step 5), so
% howards=0 there is a no-op and the check would pass vacuously. Add them when step 5 lands.
%
% beta0 MUST BE < 1 THROUGHOUT. At beta0=1, beta0beta==beta and Vhat==Vunderbar identically, so
% sourcing Vunderbar from Vhat is indistinguishable from sourcing it from itself and the bug is
% invisible. The bank's beta0=1 cross-test is structurally blind to this; that is why this file
% exists separately.

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn_nod=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
ReturnFn_d=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);

vfoptionsQH=struct();
vfoptionsQH.exoticpreferences='QuasiHyperbolic';
vfoptionsQH.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;

% The bar for the ratio check below. The two residuals it compares are both limited by the same
% cycling, so a correct Howards leaves the ratio near one (0.93 to 1.03 on the 2026-09-21 run).
% A factor of two is a judgement call, not a derived bound: there is no pre-fix measurement of
% this statistic to calibrate against, because the test was written after the fix. It is set
% wide enough that ordinary run-to-run movement cannot trip it, and tight enough that a
% Vunderbar driven toward the wrong fixed point would. Re-derive it if a pre-fix baseline is
% ever measured; do not quietly widen it to make a run pass.
ratiobar=2;

fprintf('\n===== Howards on vs off (Sophisticated) ===== \n')
fprintf('beta0=%2.2f (must be <1, or the bug being tested for is invisible) \n',Params.beta0)

%% Positive control: Naive, where Howards sits in the standard exponential leg
% Naive's Howards is the ordinary exponential-discounting one, which is known good and whose
% iteration DOES converge. So here on vs off must agree to convergence tolerance. If this line is
% not near zero, the harness itself is wrong and nothing below can be read.
vfN=vfoptionsQH; vfN.quasi_hyperbolic='Naive';
vfN_noH=vfN; vfN_noH.howards=0;
[VNon,PNon,VNonalt]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,[],vfN);
[VNoff,PNoff,VNoffalt]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,[],vfN_noH);
fprintf('CONTROL Naive nod: howards on vs off, V, this should be zero to convergence tolerance: %2.8e \n',max(abs(VNon(:)-VNoff(:))))
fprintf('CONTROL Naive nod: howards on vs off, Valt, this should be zero to convergence tolerance: %2.8e \n',max(abs(VNonalt(:)-VNoffalt(:))))
fprintf('CONTROL Naive nod: howards on vs off, Policy, this should be zero: %.3e \n',max(abs(PNon(:)-PNoff(:))))
clear VNon PNon VNonalt VNoff PNoff VNoffalt

%% (a) Sophisticated, nod, pure discretization -> S_nod_raw
fprintf('--- (a) Sophisticated, nod, pure discretization (S_nod_raw) --- \n')
vfS=vfoptionsQH; vfS.quasi_hyperbolic='Sophisticated';
vfS_noH=vfS; vfS_noH.howards=0;

[VaOn,PaOn,VuaOn]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,[],vfS);
[VaOff,PaOff,VuaOff]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,[],vfS_noH);

[~,VuaOnTrue]=ValueFnFromPolicy_InfHorz(PaOn,0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,vfS);
[~,VuaOffTrue]=ValueFnFromPolicy_InfHorz(PaOff,0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,vfS_noH);
residaOn=max(abs(VuaOn(:)-VuaOnTrue(:)));
residaOff=max(abs(VuaOff(:)-VuaOffTrue(:)));
fprintf('(a) self-consistency of Vunderbar with howards ON  (vs its own policy), diagnostic: %2.8e \n',residaOn)
fprintf('(a) self-consistency of Vunderbar with howards OFF (vs its own policy), diagnostic: %2.8e \n',residaOff)
ratioa=gather(residaOn/residaOff);
fprintf('(a) ratio ON/OFF, diagnostic (materially above one is the Howards bug): %2.4f \n',ratioa)
fprintf('(a) Howards does not bias Vunderbar relative to its own policy, this should be one: %d \n',ratioa<ratiobar)
fprintf('(a) howards on vs off, Vhat, CONVERGENCE-DEPENDENT (not a zero-check): %2.8e \n',max(abs(VaOn(:)-VaOff(:))))
fprintf('(a) howards on vs off, Policy, CONVERGENCE-DEPENDENT: %.3e \n',max(abs(PaOn(:)-PaOff(:))))
clear VaOn PaOn VuaOn VaOff PaOff VuaOff VuaOnTrue VuaOffTrue

%% (b) Sophisticated, d, refinement (the default) -> S_nod_raw
fprintf('--- (b) Sophisticated, d, refinement (S_nod_raw via ..._QuasiHyperbolic_Refine) --- \n')
[VbOn,PbOn,VubOn]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,[],vfS);
[VbOff,PbOff,VubOff]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,[],vfS_noH);

[~,VubOnTrue]=ValueFnFromPolicy_InfHorz(PbOn,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,vfS);
[~,VubOffTrue]=ValueFnFromPolicy_InfHorz(PbOff,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,vfS_noH);
residbOn=max(abs(VubOn(:)-VubOnTrue(:)));
residbOff=max(abs(VubOff(:)-VubOffTrue(:)));
fprintf('(b) self-consistency of Vunderbar with howards ON  (vs its own policy), diagnostic: %2.8e \n',residbOn)
fprintf('(b) self-consistency of Vunderbar with howards OFF (vs its own policy), diagnostic: %2.8e \n',residbOff)
ratiob=gather(residbOn/residbOff);
fprintf('(b) ratio ON/OFF, diagnostic (materially above one is the Howards bug): %2.4f \n',ratiob)
fprintf('(b) Howards does not bias Vunderbar relative to its own policy, this should be one: %d \n',ratiob<ratiobar)
fprintf('(b) howards on vs off, Vhat, CONVERGENCE-DEPENDENT (not a zero-check): %2.8e \n',max(abs(VbOn(:)-VbOff(:))))
fprintf('(b) howards on vs off, Policy, CONVERGENCE-DEPENDENT: %.3e \n',max(abs(PbOn(:)-PbOff(:))))
clear VbOn PbOn VubOn VbOff PbOff VubOff VubOnTrue VubOffTrue

%% (c) Sophisticated, d, pure discretization -> S_raw. THE PATH THAT CARRIED THE LIVE BUG.
fprintf('--- (c) Sophisticated, d, pure discretization (S_raw -- the path that carried the live bug) --- \n')
vfSc=vfS; vfSc.solnmethod='purediscretization';
vfSc_noH=vfSc; vfSc_noH.howards=0;
[VcOn,PcOn,VucOn]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,[],vfSc);
[VcOff,PcOff,VucOff]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,[],vfSc_noH);

[~,VucOnTrue]=ValueFnFromPolicy_InfHorz(PcOn,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,vfSc);
[~,VucOffTrue]=ValueFnFromPolicy_InfHorz(PcOff,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,vfSc_noH);
residcOn=max(abs(VucOn(:)-VucOnTrue(:)));
residcOff=max(abs(VucOff(:)-VucOffTrue(:)));
fprintf('(c) self-consistency of Vunderbar with howards ON  (vs its own policy), diagnostic: %2.8e \n',residcOn)
fprintf('(c) self-consistency of Vunderbar with howards OFF (vs its own policy), diagnostic: %2.8e \n',residcOff)
ratioc=gather(residcOn/residcOff);
fprintf('(c) ratio ON/OFF, diagnostic (materially above one is the Howards bug): %2.4f \n',ratioc)
fprintf('(c) Howards does not bias Vunderbar relative to its own policy, this should be one: %d \n',ratioc<ratiobar)
fprintf('(c) howards on vs off, Vhat, CONVERGENCE-DEPENDENT (not a zero-check): %2.8e \n',max(abs(VcOn(:)-VcOff(:))))
fprintf('(c) howards on vs off, Policy, CONVERGENCE-DEPENDENT: %.3e \n',max(abs(PcOn(:)-PcOff(:))))
clear VcOn PcOn VucOn VcOff PcOff VucOff VucOnTrue VucOffTrue

%% (d) and (e): the two postGI Sophisticated raws, which got Howards on 2026-10-21 (step 5)
% These are NEW and UNVALIDATED solver code, and they are the hard case: the interpolated
% continuation at a fixed fine index is a two-point weighted combination of coarse EV values, so
% the block resolves the window position by hand rather than reusing the raw's interp1. Getting
% that index arithmetic wrong is the obvious failure mode and these two checks are what would
% catch it -- without them step 5 ships untested.
% Same self-consistency statistic as (a)-(c): ValueFnFromPolicy re-solves the Vunderbar fixed point
% independently, so a Howards block that drives Vunderbar somewhere it should not go shows up as a
% worse ON residual than OFF. Note these two paths are ALSO limited by the kinked-Vunderbar
% collapse at beta0<1, which Howards does not touch and which affects ON and OFF equally -- so it
% cancels in the ratio and does not interfere with what is being tested here.
for giarm=1:2
    if giarm==1
        armlabel='(d) Sophisticated, nod, grid interpolation (S_postGI_nod_raw)';
        n_darm=0; d_gridarm=[]; ReturnFnarm=ReturnFn_nod;
    else
        armlabel='(e) Sophisticated, d, refine + grid interpolation (S_Refine_postGI_raw)';
        n_darm=n_d; d_gridarm=d_grid; ReturnFnarm=ReturnFn_d;
    end
    vfG=vfS; vfG.gridinterplayer=1; vfG.ngridinterp=5; vfG.howardssparse=0;
    % No default for Sophisticated + GI; set it to the recommended value.
    vfG.maxaprimediff=ceil(N_a/(Params.beta0*4));
    vfG_noH=vfG; vfG_noH.howards=0;

    [VgOn,PgOn,VugOn]=ValueFnIter_InfHorz(n_darm,n_a,n_z,d_gridarm,a_grid,z_grid,pi_z,ReturnFnarm,Params,DiscountFactorParamNames,[],vfG);
    [VgOff,PgOff,VugOff]=ValueFnIter_InfHorz(n_darm,n_a,n_z,d_gridarm,a_grid,z_grid,pi_z,ReturnFnarm,Params,DiscountFactorParamNames,[],vfG_noH);

    [~,VugOnTrue]=ValueFnFromPolicy_InfHorz(PgOn,n_darm,n_a,n_z,d_gridarm,a_grid,z_grid,pi_z,ReturnFnarm,Params,DiscountFactorParamNames,vfG);
    [~,VugOffTrue]=ValueFnFromPolicy_InfHorz(PgOff,n_darm,n_a,n_z,d_gridarm,a_grid,z_grid,pi_z,ReturnFnarm,Params,DiscountFactorParamNames,vfG_noH);
    residgOn=max(abs(VugOn(:)-VugOnTrue(:)));
    residgOff=max(abs(VugOff(:)-VugOffTrue(:)));
    ratiog=gather(residgOn/residgOff);

    fprintf('--- %s --- \n',armlabel)
    fprintf('    self-consistency of Vunderbar with howards ON  (vs its own policy), diagnostic: %2.8e \n',residgOn)
    fprintf('    self-consistency of Vunderbar with howards OFF (vs its own policy), diagnostic: %2.8e \n',residgOff)
    fprintf('    ratio ON/OFF, diagnostic (materially above one means the new Howards block is wrong): %2.4f \n',ratiog)
    fprintf('    postGI Howards does not bias Vunderbar relative to its own policy, this should be one: %d \n',ratiog<ratiobar)
    fprintf('    howards on vs off, Vhat, CONVERGENCE-DEPENDENT (not a zero-check): %2.8e \n',max(abs(VgOn(:)-VgOff(:))))
    clear VgOn PgOn VugOn VgOff PgOff VugOff VugOnTrue VugOffTrue
end

%% Near-exponential control: beta0 close to 1, where the iteration should actually converge
% At beta0=0.99 the quasi-hyperbolic problem is a small perturbation of the exponential one, so the
% Sophisticated policy should settle rather than cycle. If it does converge (no maxiter warning
% from the raw, and the self-consistency residuals below are near zero), then howards on vs off
% becomes a genuine EXACT check, which is what a plain accelerator test wants.
% beta0beta=0.99*beta is still clearly distinct from beta, so the old bug -- whose fixed point
% carries beta*beta0beta and an extra flow payoff -- would remain plainly visible here.
% If the residuals below are NOT near zero, this case did not converge either; read the on-vs-off
% lines as convergence-dependent and treat the ratio checks above as the real test.
fprintf('--- near-exponential control, beta0=0.99, d + pure discretization (S_raw) --- \n')
Params_nearexp=Params;
Params_nearexp.beta0=0.99;
[VdOn,PdOn,VudOn]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params_nearexp,DiscountFactorParamNames,[],vfSc);
[VdOff,PdOff,VudOff]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params_nearexp,DiscountFactorParamNames,[],vfSc_noH);

[~,VudOnTrue]=ValueFnFromPolicy_InfHorz(PdOn,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params_nearexp,DiscountFactorParamNames,vfSc);

% NOTE ON THE CHECK COUNT. This block emits 3 identity checks when the gate passes and 0 when it
% does not (the diagnostics that replace them carry no 'should be' wording, so CoreSummary ignores
% them). The bank's total check count therefore moves by 3 depending on model behaviour rather
% than on code. That is deliberate -- the alternative is to assert an identity that is not true --
% but it means a 3-check drop here is NOT evidence that a check stopped being executed. The
% CONVERGED / did NOT converge line printed below says which branch ran; read it before concluding
% anything from the count.
%
% Whether on-vs-off is an IDENTITY here depends on whether this case converged, so the check is
% gated on that rather than asserted and hedged in prose. Claiming 'should be zero' on a run that
% stopped at maxiter would be claiming something untrue, and a diary summary reads the wording, not
% the hedge. Converged means the returned Vunderbar IS the policy-evaluation value of its own
% policy; with tolerance 1e-9 and the 1/(1-beta) amplification that is about 2e-8, so 1e-6 sits
% clear of it and four orders below the 1.4e-02 this case produced on 2026-09-21.
convresid=gather(max(abs(VudOn(:)-VudOnTrue(:))));
fprintf('beta0=0.99: self-consistency of Vunderbar with howards ON, the convergence gate: %2.8e \n',convresid)
if convresid<10^(-6)
    fprintf('beta0=0.99: CONVERGED, so howards on vs off is a genuine identity and is claimed as one \n')
    fprintf('beta0=0.99: howards on vs off, Vhat, this should be zero: %2.8e \n',max(abs(VdOn(:)-VdOff(:))))
    fprintf('beta0=0.99: howards on vs off, Vunderbar, this should be zero: %2.8e \n',max(abs(VudOn(:)-VudOff(:))))
    fprintf('beta0=0.99: howards on vs off, Policy, this should be zero: %.3e \n',max(abs(PdOn(:)-PdOff(:))))
else
    fprintf('beta0=0.99: did NOT converge, so howards on vs off is NOT an identity here \n')
    fprintf('beta0=0.99: howards on vs off, Vhat, CONVERGENCE-DEPENDENT: %2.8e \n',max(abs(VdOn(:)-VdOff(:))))
    fprintf('beta0=0.99: howards on vs off, Vunderbar, CONVERGENCE-DEPENDENT: %2.8e \n',max(abs(VudOn(:)-VudOff(:))))
    fprintf('beta0=0.99: howards on vs off, Policy, CONVERGENCE-DEPENDENT: %.3e \n',max(abs(PdOn(:)-PdOff(:))))
end

%% Anti-vacuity: beta0 must actually differ from 1 in both cases used above
fprintf('anti-vacuity: beta0=%2.2f and beta0=0.99 both give beta0*beta distinct from beta (%2.6f, %2.6f vs %2.6f) \n', ...
    Params.beta0,Params.beta0*Params.beta,0.99*Params.beta,Params.beta)

%%
output=struct();

end
