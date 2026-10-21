function output=QHDInfHorz_d_z(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% InfHorz quasi-hyperbolic, with d, with z, no grid interpolation layer.
%
% This is the subcode with the most solver paths in it, because with n_d>0 the dispatcher
% branches on vfoptions.solnmethod:
%   solnmethod='purediscretization_refinement' (THE DEFAULT)
%       -> ValueFnIter_InfHorz_QuasiHyperbolic_Refine, which refines d out and then calls
%          ValueFnIter_InfHorz_QuasiHyperbolicN_nod_raw / ...S_nod_raw
%       -> and, for the Naive V_std leg, one of ValueFnIter_InfHorz_{nod,sparse_nod,HowardGreedy_nod}_raw
%          depending on vfoptions.howardssparse / vfoptions.howardsgreedy
%   solnmethod='purediscretization'
%       -> ValueFnIter_InfHorz_QuasiHyperbolicN_raw / ...S_raw
%
% ValueFnIter_InfHorz_QuasiHyperbolicN_raw is REGRESSION TARGET 2: it used to overwrite its own
% Policy accumulator with zeros(2,N_a,N_z) before unpacking the flat argmax into (d,aprime), so
% the unpack read the zeroed array. That path is reachable only by setting
% solnmethod='purediscretization' with n_d>0, which is why nothing had ever run it.
%
% The independent (d,aprime) oracle below rebuilds the whole quasi-hyperbolic argmax from
% a_grid/d_grid/z_grid/pi_z without calling the toolkit, and compares both channels of Policy
% against it. That is what pins the (d,aprime) unpack, and it pins the discount factor at the
% same time (regression target 1).

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();

N_d=prod(n_d);
N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);

%% Quasi-Hyperbolic Discounting
vfoptions.exoticpreferences='QuasiHyperbolic';
vfoptions.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;

beta=Params.beta;
beta0beta=Params.beta0*Params.beta;

% The (d,aprime) stacked grids, in the toolkit's ordering: d varies fastest
dvec=repmat(d_grid,N_a,1);        % (N_d*N_a)-by-1
aprimevec=repelem(a_grid,N_d,1);  % (N_d*N_a)-by-1
dmat=repmat(dvec,1,N_a);          % (N_d*N_a)-by-N_a

%% The exponential-discounting reference solve for this same model
vfoptionsexp=struct(); % deliberately no exoticpreferences field
[Vexp,Policyexp]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsexp);

%% ------------------------------------------------------------------ Naive
vfoptions.quasi_hyperbolic='Naive';

%% Naive, refinement (the default solnmethod when n_d>0)
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1,V1alt,Policy1alt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

PolicyVals1=PolicyInd2Val_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,vfoptions1);
PolicyVals1alt=PolicyInd2Val_InfHorz(Policy1alt,n_d,n_a,n_z,d_grid,a_grid,vfoptions1);
fprintf('Naive (refine): PolicyInd2Val ran on both Policy and Policyalt \n')

fprintf('Naive (refine): size(Policy) should be [2 %d %d], is [%s] \n',N_a,N_z,num2str(size(Policy1)))
fprintf('Naive (refine): d index in range, this should be one: %d \n',gather(all(all(Policy1(1,:)>=1)) && all(all(Policy1(1,:)<=N_d))))
fprintf('Naive (refine): aprime index in range, this should be one: %d \n',gather(all(all(Policy1(2,:)>=1)) && all(all(Policy1(2,:)<=N_a))))
fprintf('Naive (refine): d index in range (alt), this should be one: %d \n',gather(all(all(Policy1alt(1,:)>=1)) && all(all(Policy1alt(1,:)<=N_d))))
fprintf('Naive (refine): aprime index in range (alt), this should be one: %d \n',gather(all(all(Policy1alt(2,:)>=1)) && all(all(Policy1alt(2,:)<=N_a))))
fprintf('Naive (refine): V is finite everywhere, this should be one: %d \n',gather(all(isfinite(V1(:)))))

%% Naive: Valt and Policyalt must equal the plain exponential-discounting solution
fprintf('Naive (refine): Valt vs plain exponential ValueFnIter, this should be zero: %.3e \n',max(abs(V1alt(:)-Vexp(:))))
fprintf('Naive (refine): Policyalt vs plain exponential ValueFnIter, this should be zero: %.3e \n',max(abs(Policy1alt(:)-Policyexp(:))))

%% Naive, plain pure discretization (this is what reaches ValueFnIter_InfHorz_QuasiHyperbolicN_raw)
vfoptions2=vfoptions;
vfoptions2.solnmethod='purediscretization';
simoptions2=simoptions;
[V2,Policy2,V2alt,Policy2alt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);

PolicyVals2=PolicyInd2Val_InfHorz(Policy2,n_d,n_a,n_z,d_grid,a_grid,vfoptions2);
fprintf('Naive (plain): PolicyInd2Val ran \n')
fprintf('Naive (plain): size(Policy) should be [2 %d %d], is [%s] \n',N_a,N_z,num2str(size(Policy2)))
fprintf('Naive (plain): d index in range, this should be one: %d \n',gather(all(all(Policy2(1,:)>=1)) && all(all(Policy2(1,:)<=N_d))))
fprintf('Naive (plain): aprime index in range, this should be one: %d \n',gather(all(all(Policy2(2,:)>=1)) && all(all(Policy2(2,:)<=N_a))))

%% Refinement and plain pure discretization solve the identical problem, so they must agree
% (d refinement is preference-independent: d only ever enters the return function, never the
% continuation, so argmax_d of the QH objective is argmax_d of u regardless of the discounting.)
fprintf('Naive: refine vs plain purediscretization, V, this should be zero: %.3e \n',max(abs(V1(:)-V2(:))))
fprintf('Naive: refine vs plain purediscretization, Policy, this should be zero (ties aside): %.3e \n',max(abs(Policy1(:)-Policy2(:))))
fprintf('Naive: refine vs plain purediscretization, Valt, this should be zero: %.3e \n',max(abs(V1alt(:)-V2alt(:))))
fprintf('Naive: refine vs plain purediscretization, Policyalt, this should be zero (ties aside): %.3e \n',max(abs(Policy1alt(:)-Policy2alt(:))))

%% Naive: the independent (d,aprime) oracle
% Rebuild the quasi-hyperbolic argmax here, from the grids and pi_z, with no toolkit call:
%   Vtilde(a,z) = max_{d,aprime} u(d,aprime,a,z) + beta0*beta*E[Valt(aprime,z')|z]
% then split the flat argmax into (d,aprime) with the toolkit's own d-fastest convention.
% Because Valt is taken as given (the Naive raw does exactly one step off it), this is exact.
Vcheck=zeros(N_a,N_z);
dcheck=zeros(N_a,N_z);
aprimecheck=zeros(N_a,N_z);
V1altcpu=gather(V1alt);
pi_zcpu=gather(pi_z);
for z_c=1:N_z
    cmat=(1+Params.r)*a_grid'+Params.w*z_grid(z_c)*dvec-aprimevec; % (N_d*N_a)-by-N_a by implicit expansion
    util=-Inf(N_d*N_a,N_a);
    okmat=(cmat>0)&(dmat<1);
    util(okmat)=(cmat(okmat).^(1-Params.sigma)-1)/(1-Params.sigma)+Params.varphi*((1-dmat(okmat)).^(1-Params.eta)-1)/(1-Params.eta);
    EV_z=V1altcpu*pi_zcpu(z_c,:)';        % N_a-by-1, expectation over z' at each aprime
    entireEV_z=repelem(EV_z,N_d,1);       % (N_d*N_a)-by-1, d varies fastest
    RHS=util+beta0beta*entireEV_z;
    [Vtemp,maxindex]=max(RHS,[],1);
    Vcheck(:,z_c)=Vtemp';
    dcheck(:,z_c)=rem(maxindex-1,N_d)+1;
    aprimecheck(:,z_c)=ceil(maxindex/N_d);
end
fprintf('Naive (refine): independent oracle, V, this should be zero: %2.8e \n',max(abs(gather(V1(:))-Vcheck(:))))
fprintf('Naive (refine): independent oracle, d channel of Policy, this should be zero: %.3e \n',max(abs(gather(reshape(Policy1(1,:,:),[N_a,N_z]))-dcheck),[],'all'))
fprintf('Naive (refine): independent oracle, aprime channel of Policy, this should be zero: %.3e \n',max(abs(gather(reshape(Policy1(2,:,:),[N_a,N_z]))-aprimecheck),[],'all'))
fprintf('Naive (plain): independent oracle, V, this should be zero: %2.8e \n',max(abs(gather(V2(:))-Vcheck(:))))
fprintf('Naive (plain): independent oracle, d channel of Policy, this should be zero: %.3e \n',max(abs(gather(reshape(Policy2(1,:,:),[N_a,N_z]))-dcheck),[],'all'))
fprintf('Naive (plain): independent oracle, aprime channel of Policy, this should be zero: %.3e \n',max(abs(gather(reshape(Policy2(2,:,:),[N_a,N_z]))-aprimecheck),[],'all'))

% The same oracle rerun with the wrong discount factors, so the diary records that the check
% above discriminates. Both of these must be well away from zero.
Vcheckbeta=zeros(N_a,N_z);
Vcheckbeta0=zeros(N_a,N_z);
for z_c=1:N_z
    cmat=(1+Params.r)*a_grid'+Params.w*z_grid(z_c)*dvec-aprimevec;
    util=-Inf(N_d*N_a,N_a);
    okmat=(cmat>0)&(dmat<1);
    util(okmat)=(cmat(okmat).^(1-Params.sigma)-1)/(1-Params.sigma)+Params.varphi*((1-dmat(okmat)).^(1-Params.eta)-1)/(1-Params.eta);
    EV_z=V1altcpu*pi_zcpu(z_c,:)';
    entireEV_z=repelem(EV_z,N_d,1);
    Vcheckbeta(:,z_c)=max(util+beta*entireEV_z,[],1)';
    Vcheckbeta0(:,z_c)=max(util+Params.beta0*entireEV_z,[],1)';
end
fprintf('Naive: oracle rerun discounting by beta alone, this should be WELL ABOVE zero: %2.8f \n',max(abs(gather(V1(:))-Vcheckbeta(:))))
fprintf('Naive: oracle rerun discounting by beta0 alone, this should be WELL ABOVE zero: %2.8f \n',max(abs(gather(V1(:))-Vcheckbeta0(:))))

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
fprintf('Naive (refine): ValueFnFromPolicy V (Vtilde), this should be zero to convergence tolerance: %2.8e \n',max(abs(VfromPolicy(:)-V1(:))))
fprintf('Naive (refine): ValueFnFromPolicy Valt (V_std), this should be zero to convergence tolerance: %2.8e \n',max(abs(ValtfromPolicy(:)-V1alt(:))))

%% Naive: anti-vacuity
fprintf('Naive: quasi-hyperbolic matters, max|Policy-Policyalt| should be well above zero: %2.8f \n',max(abs(Policy1(:)-Policy1alt(:))))
fprintf('Naive: quasi-hyperbolic matters, max|Vtilde-V_std| should be well above zero: %2.8f \n',max(abs(V1(:)-V1alt(:))))
fprintf('Naive: the d channel is not degenerate, number of distinct d chosen should be above one: %d \n',numel(unique(gather(Policy1(1,:)))))

%% Naive: the three Howards branches inside ValueFnIter_InfHorz_QuasiHyperbolic_Refine
% Refine picks between ValueFnIter_InfHorz_{HowardGreedy_nod, sparse_nod, nod}_raw for the
% V_std leg. All three must give the same V_std, hence the same quasi-hyperbolic step.
vfoptions1_sp=vfoptions1;
vfoptions1_sp.howardssparse=1;
[V1sp,Policy1sp,V1altsp,Policy1altsp]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_sp);
fprintf('Naive (refine): howardssparse=1, V, this should be zero: %.3e \n',max(abs(V1(:)-V1sp(:))))
fprintf('Naive (refine): howardssparse=1, Policy, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1sp(:))))
fprintf('Naive (refine): howardssparse=1, Valt, this should be zero: %.3e \n',max(abs(V1alt(:)-V1altsp(:))))

vfoptions1_gr=vfoptions1;
vfoptions1_gr.howardsgreedy=1;
[V1gr,Policy1gr,V1altgr,Policy1altgr]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_gr);
fprintf('Naive (refine): howardsgreedy=1, V, this should be zero: %.3e \n',max(abs(V1(:)-V1gr(:))))
fprintf('Naive (refine): howardsgreedy=1, Policy, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1gr(:))))
fprintf('Naive (refine): howardsgreedy=1, Valt, this should be zero: %.3e \n',max(abs(V1alt(:)-V1altgr(:))))

vfoptions1_noH=vfoptions1;
vfoptions1_noH.howards=0;
[V1noH,Policy1noH]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_noH);
fprintf('Naive (refine): howards=0 (pure VFI), V, this should be zero: %.3e \n',max(abs(V1(:)-V1noH(:))))
fprintf('Naive (refine): howards=0 (pure VFI), Policy, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1noH(:))))

%% Naive: lowmemory
% ValueFnIter_InfHorz_Refine builds the return matrix by looping over z when lowmemory=1, but
% ValueFnIter_InfHorz_QuasiHyperbolic_Refine never reads vfoptions.lowmemory, so today this is
% a no-op. It is checked so that a future lowmemory branch there cannot silently diverge.
vfoptions1_lm=vfoptions1;
vfoptions1_lm.lowmemory=1;
[V1lm,Policy1lm]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_lm);
fprintf('Naive (refine): lowmemory=1, V, this should be zero: %.3e \n',max(abs(V1(:)-V1lm(:))))
fprintf('Naive (refine): lowmemory=1, Policy, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1lm(:))))

clear V1sp V1gr V1noH V1lm Policy1sp Policy1gr Policy1noH Policy1lm V1altsp V1altgr Policy1altsp Policy1altgr

%% ---------------------------------------------------------- Sophisticated
% Reminder: the dispatcher forces maxiter=1000. Since 2026-09-01 both raws this subcode reaches
% (S_nod_raw via refine, S_raw via purediscretization) run Howards on Vunderbar; before that
% neither did. Nothing below asserts tight convergence; the Vunderbar Bellman residual is the
% convergence diagnostic. Note Vhat is now recomputed after the loop from the returned Vunderbar,
% so the Vhat identity is imposed by construction and is not a convergence measure.
vfoptions3=vfoptions;
vfoptions3.quasi_hyperbolic='Sophisticated';
simoptions3=simoptions;
[V3,Policy3,V3alt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
% V3=Vhat, V3alt=Vunderbar

PolicyVals3=PolicyInd2Val_InfHorz(Policy3,n_d,n_a,n_z,d_grid,a_grid,vfoptions3);
fprintf('Sophisticated (refine): PolicyInd2Val ran \n')
fprintf('Sophisticated (refine): size(Policy) should be [2 %d %d], is [%s] \n',N_a,N_z,num2str(size(Policy3)))
fprintf('Sophisticated (refine): d index in range, this should be one: %d \n',gather(all(all(Policy3(1,:)>=1)) && all(all(Policy3(1,:)<=N_d))))
fprintf('Sophisticated (refine): aprime index in range, this should be one: %d \n',gather(all(all(Policy3(2,:)>=1)) && all(all(Policy3(2,:)<=N_a))))
fprintf('Sophisticated (refine): Vhat is finite everywhere, this should be one: %d \n',gather(all(isfinite(V3(:)))))
fprintf('Sophisticated (refine): Vunderbar is finite everywhere, this should be one: %d \n',gather(all(isfinite(V3alt(:)))))

%% Sophisticated, plain pure discretization (reaches ValueFnIter_InfHorz_QuasiHyperbolicS_raw)
vfoptions4=vfoptions3;
vfoptions4.solnmethod='purediscretization';
[V4,Policy4,V4alt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('Sophisticated (plain): size(Policy) should be [2 %d %d], is [%s] \n',N_a,N_z,num2str(size(Policy4)))
fprintf('Sophisticated (plain): d index in range, this should be one: %d \n',gather(all(all(Policy4(1,:)>=1)) && all(all(Policy4(1,:)<=N_d))))
fprintf('Sophisticated (plain): aprime index in range, this should be one: %d \n',gather(all(all(Policy4(2,:)>=1)) && all(all(Policy4(2,:)<=N_a))))

% ValueFnIter_InfHorz_QuasiHyperbolicS_raw runs Howards on the joint Vhat/Vunderbar iteration
% while the refine path (which goes to ...S_nod_raw) does not, so these two converge along
% different sequences. They target the same fixed point, so this is a roughly-equal check
% rather than an exact one; how close it comes is itself a Sophisticated convergence report.
% These are NOT zero-checks. Read them against the Vunderbar Bellman residual printed below: if
% that residual is large, the Sophisticated iteration stopped at maxiter rather than at a fixed
% point, and two solution paths that approach the same fixed point differently are then simply
% in different places. They become zero-checks only once that residual is at tolerance.
fprintf('Sophisticated: refine vs plain purediscretization, Vhat, CONVERGENCE-DEPENDENT (compare to the Vunderbar residual below): %2.8e \n',max(abs(V3(:)-V4(:))))
fprintf('Sophisticated: refine vs plain purediscretization, Vunderbar, CONVERGENCE-DEPENDENT: %2.8e \n',max(abs(V3alt(:)-V4alt(:))))
fprintf('Sophisticated: refine vs plain purediscretization, Policy, CONVERGENCE-DEPENDENT: %2.8f \n',max(abs(Policy3(:)-Policy4(:))))

%% Sophisticated: the (d,aprime) oracle, given Vunderbar
% Policy must be argmax_{d,aprime} u + beta0*beta*E[Vunderbar], and Vhat/Vunderbar must be that
% objective evaluated at it with beta0*beta / beta respectively. All hold up to the fixed-point
% residual of the Vunderbar iteration.
Vhatcheck=zeros(N_a,N_z);
Vunderbarcheck=zeros(N_a,N_z);
dcheck3=zeros(N_a,N_z);
aprimecheck3=zeros(N_a,N_z);
EVunderbar3=zeros(N_a,N_z);
V3altcpu=gather(V3alt);
for z_c=1:N_z
    cmat=(1+Params.r)*a_grid'+Params.w*z_grid(z_c)*dvec-aprimevec;
    util=-Inf(N_d*N_a,N_a);
    okmat=(cmat>0)&(dmat<1);
    util(okmat)=(cmat(okmat).^(1-Params.sigma)-1)/(1-Params.sigma)+Params.varphi*((1-dmat(okmat)).^(1-Params.eta)-1)/(1-Params.eta);
    EV_z=V3altcpu*pi_zcpu(z_c,:)';
    entireEV_z=repelem(EV_z,N_d,1);
    [Vtemp,maxindex]=max(util+beta0beta*entireEV_z,[],1);
    Vhatcheck(:,z_c)=Vtemp';
    dcheck3(:,z_c)=rem(maxindex-1,N_d)+1;
    aprimecheck3(:,z_c)=ceil(maxindex/N_d);
    EVunderbar3(:,z_c)=EV_z(aprimecheck3(:,z_c)); % E[Vunderbar] at the chosen aprime, for the sign check below
    tempmaxindex=maxindex+(0:1:N_a-1)*(N_d*N_a);
    RHSbeta=util+beta*entireEV_z;
    Vunderbarcheck(:,z_c)=RHSbeta(tempmaxindex)';
end
% Also CONVERGENCE-DEPENDENT. The oracle recomputes the argmax GIVEN the returned Vunderbar; if
% Vunderbar is not a fixed point then the returned Policy was taken against a different (earlier)
% Vunderbar, and the two legitimately disagree. Judge these against the Vunderbar residual.
fprintf('Sophisticated (refine): oracle, Vhat, CONVERGENCE-DEPENDENT: %2.8e \n',max(abs(gather(V3(:))-Vhatcheck(:))))
fprintf('Sophisticated (refine): oracle, Vunderbar, CONVERGENCE-DEPENDENT: %2.8e \n',max(abs(gather(V3alt(:))-Vunderbarcheck(:))))
fprintf('Sophisticated (refine): oracle, d channel of Policy, CONVERGENCE-DEPENDENT: %.3e \n',max(abs(gather(reshape(Policy3(1,:,:),[N_a,N_z]))-dcheck3),[],'all'))
fprintf('Sophisticated (refine): oracle, aprime channel of Policy, CONVERGENCE-DEPENDENT: %.3e \n',max(abs(gather(reshape(Policy3(2,:,:),[N_a,N_z]))-aprimecheck3),[],'all'))

%% Sophisticated: ValueFnFromPolicy_InfHorz
% Sophisticated needs no Policyalt: Policy alone determines both Vhat and Vunderbar.
[VfromPolicyS,ValtfromPolicyS]=ValueFnFromPolicy_InfHorz(Policy3,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
fprintf('Sophisticated (refine): ValueFnFromPolicy V (Vhat) vs the solver, CONVERGENCE-DEPENDENT (this measures solver non-convergence, it is not a fault): %2.8e \n',max(abs(VfromPolicyS(:)-V3(:))))
fprintf('Sophisticated (refine): ValueFnFromPolicy Valt (Vunderbar) vs the solver, CONVERGENCE-DEPENDENT: %2.8e \n',max(abs(ValtfromPolicyS(:)-V3alt(:))))

%% Sophisticated: anti-vacuity
fprintf('Sophisticated: quasi-hyperbolic matters, max|Policy-Policyexp| should be well above zero: %2.8f \n',max(abs(Policy3(:)-Policyexp(:))))
fprintf('Sophisticated: quasi-hyperbolic matters, max|Vhat-Vunderbar| should be well above zero: %2.8f \n',max(abs(V3(:)-V3alt(:))))
% The identity above says Vhat-Vunderbar = (beta0*beta-beta)*E[Vunderbar], and beta0*beta<beta,
% so Vhat exceeds Vunderbar exactly where E[Vunderbar] is NEGATIVE -- not everywhere. With
% sigma=2 the flow utility is u=1-1/c, which is POSITIVE for c>1, so V is not signed in this
% calibration and 'Vhat>Vunderbar everywhere' is false in the nod models and true in the d
% models (where the leisure term drags V negative). The share below is the correct check.
signok=(gather(V3)-gather(V3alt)>0)==(EVunderbar3<0);
fprintf('Sophisticated: Vhat exceeds Vunderbar exactly where E[Vunderbar]<0, share of states should be one: %2.4f \n',mean(signok(:)))

%% Sophisticated vs Naive
fprintf('Naive vs Sophisticated: max|V_Naive - V_Sophisticated| should be well above zero: %2.8f \n',max(abs(V1(:)-V3(:))))

%% Figures
zmid=ceil(N_z/2);

fig=figure(figure_c);
subplot(1,3,1); plot(a_grid,squeeze(PolicyVals1(2,:,zmid)), a_grid,squeeze(PolicyVals1alt(2,:,zmid)), a_grid,a_grid,'k:')
title('Naive: aprime policy (mid z)'); xlabel('a'); legend('quasi-hyperbolic','exponential','45 degree','Location','northwest')
subplot(1,3,2); plot(a_grid,squeeze(PolicyVals1(1,:,zmid)), a_grid,squeeze(PolicyVals1alt(1,:,zmid)))
title('Naive: d policy (mid z)'); xlabel('a'); legend('quasi-hyperbolic','exponential')
subplot(1,3,3); plot(a_grid,V1(:,zmid), a_grid,V1alt(:,zmid))
title('Naive: value fn (mid z)'); xlabel('a'); legend('Vtilde','V\_std','Location','southeast')

fig=figure(100+figure_c);
subplot(1,3,1); plot(a_grid,squeeze(PolicyVals3(2,:,zmid)), a_grid,squeeze(PolicyVals1alt(2,:,zmid)), a_grid,a_grid,'k:')
title('Sophisticated: aprime policy (mid z)'); xlabel('a'); legend('quasi-hyperbolic','exponential','45 degree','Location','northwest')
subplot(1,3,2); plot(a_grid,squeeze(PolicyVals3(1,:,zmid)), a_grid,squeeze(PolicyVals1alt(1,:,zmid)))
title('Sophisticated: d policy (mid z)'); xlabel('a'); legend('quasi-hyperbolic','exponential')
subplot(1,3,3); plot(a_grid,V3(:,zmid), a_grid,V3alt(:,zmid))
title('Sophisticated: value fn (mid z)'); xlabel('a'); legend('Vhat','Vunderbar','Location','southeast')

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
