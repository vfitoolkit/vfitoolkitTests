function output=EZInfHorz_CrossTests_collapse(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames)
% Cross-tests: the special Epstein-Zin identities. These play the role the beta0=1 cross-test
% plays in the QH bank: each is a parameter value at which Epstein-Zin must reproduce a solve
% that does not go through the Epstein-Zin raws at all.
%
% (i)   cons-units, gamma=1/phi: Epstein-Zin collapses to vNM expected utility with period utility
%       x^(1-gamma)/(1-gamma) (= EZReturnFn_negativeUtils at ezsigma=ezgamma), with
%       V_EZ = ((1-gamma)*V_vNM)^(1/(1-gamma)) and the same policy up to ties.
% (ii)  utility-units, EZriskaversion=0: Epstein-Zin collapses to vNM with the same ReturnFn,
%       positive and negative cases (tests the ezc3/ezc4 sign handling and the 1+/-crisk branch).
% (iv)  vfoptions.EZoneminusbeta=1 versus manually scaling the return fn:
%       cons-units is homogeneous of degree 1 in the return, so a (1-beta) on F^(1-1/phi) is the
%       same as scaling F by (1-beta)^(1/(1-1/phi)); utility-units V=(1-beta)*u+beta*CE is the
%       default recursion with u scaled by (1-beta). All three cases.
% (xii) homogeneity (units invariance) in cons-units: scale a_grid and w by lambda and V scales
%       by lambda with the same policy. At the baseline gamma and at gamma=20 (overflow stress:
%       the certainty-equivalent raises to the powers 1-gamma and 1/(1-gamma)).
% (i) and (ii) again with d (the with-d raw is a different file).
%
% POLICY COMPARISONS ARE NOT IDENTITIES HERE. In (i), (ii) and (iv) the two sides are two
% different fixed-point iterations of the same problem, which can legitimately pick different
% grid points where the objective is within tolerance. So the raw index difference is printed as
% a DIAGNOSTIC, and what is graded is the Epstein-Zin VALUE of the other solve's policy
% (ValueFnFromPolicy_InfHorz under the EZ vfoptions, at the vNM/scaled solve's Policy) against
% V_EZ. If the policies differ only at ties, that agrees to convergence tolerance; if either
% policy is wrong, V_EZ(other policy) is strictly below V_EZ somewhere (or above it, if V_EZ
% itself is not the max), and it does not. (xii) is the identical computation up to scale, so
% there the raw indices are graded.
%
% Scale factors enter the scaled return fns as PARAMETERS (Params.ezscalefactor...), because
% GPU arrayfun does not support anonymous functions that capture workspace variables.

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn_cons=@(aprime,a,z,r,w) EZReturnFn_cons_nod_z(aprime,a,z,r,w);
ReturnFn_posU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_positiveUtils_nod_z(aprime,a,z,r,w,ezsigma);
ReturnFn_negU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_negativeUtils_nod_z(aprime,a,z,r,w,ezsigma);
ReturnFn_cons_d=@(d,aprime,a,z,r,w,varphi) EZReturnFn_cons_d_z(d,aprime,a,z,r,w,varphi);
ReturnFn_posU_d=@(d,aprime,a,z,r,w,varphi,ezsigma) EZReturnFn_positiveUtils_d_z(d,aprime,a,z,r,w,varphi,ezsigma);
ReturnFn_negU_d=@(d,aprime,a,z,r,w,varphi,ezsigma) EZReturnFn_negativeUtils_d_z(d,aprime,a,z,r,w,varphi,ezsigma);

vftol=10^(-10); % see the driver header

vfN=struct(); % vNM: deliberately no exoticpreferences field
vfN.tolerance=vftol;

fprintf('\n===== Epstein-Zin collapse cross-tests ===== \n')

%% (i) cons-units, gamma=1/phi collapses Epstein-Zin to vNM expected utility
Pi=Params;
Pi.ezphi=1/Params.ezgamma;
Pi.ezsigma=Params.ezgamma; % the vNM reference's curvature
vfEZ1=struct();
vfEZ1.exoticpreferences='EpsteinZin';
vfEZ1.EZutils=0;
vfEZ1.EZriskaversion='ezgamma';
vfEZ1.EZeis='ezphi';
vfEZ1.tolerance=vftol;
[V1a,Policy1a]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Pi,DiscountFactorParamNames,[],vfEZ1);
[V1b,Policy1b]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_negU,Pi,DiscountFactorParamNames,[],vfN);
V1bt=((1-Pi.ezgamma)*V1b).^(1/(1-Pi.ezgamma));
fprintf('(i) nod: gamma=1/phi collapse, V after transform (relative), this should be zero to convergence tolerance: %2.8e \n',max(abs(V1a(:)-V1bt(:)))/max(abs(V1a(:))))
fprintf('(i) nod: gamma=1/phi collapse, raw Policy index difference, DIAGNOSTIC (near-ties can flip a grid point): %.3e at %i of %i states \n',max(abs(Policy1a(:)-Policy1b(:))),gather(sum(Policy1a(:)~=Policy1b(:))),N_a*N_z)
V1alt=ValueFnFromPolicy_InfHorz(Policy1b,0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Pi,DiscountFactorParamNames,vfEZ1);
fprintf('(i) nod: gamma=1/phi collapse, EZ value of the vNM policy vs V_EZ (relative), this should be zero to convergence tolerance: %2.8e \n',max(abs(V1a(:)-V1alt(:)))/max(abs(V1a(:))))
% and on the lowmemory path (its raw converts z to gridvals itself)
vfEZ1lm=vfEZ1;
vfEZ1lm.lowmemory=1;
[V1c,Policy1c]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Pi,DiscountFactorParamNames,[],vfEZ1lm);
fprintf('(i) nod, lowmemory=1: gamma=1/phi collapse, V after transform (relative), this should be zero to convergence tolerance: %2.8e \n',max(abs(V1c(:)-V1bt(:)))/max(abs(V1c(:))))
fprintf('(i) nod, lowmemory=1: vs lowmemory=0, Policy, this should be zero: %.3e \n',max(abs(Policy1c(:)-Policy1a(:))))
% anti-vacuity: away from gamma=1/phi the transform must NOT hold (so the line above is a test)
[V1d,~]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfEZ1); % phi=0.5, not 1/gamma
fprintf('(i) nod: the transform at phi=%g (not 1/gamma=%g), relative, should be WELL ABOVE zero: %2.8f \n',Params.ezphi,1/Params.ezgamma,max(abs(V1d(:)-V1bt(:)))/max(abs(V1d(:))))
clear V1a V1b V1bt V1c V1d V1alt Policy1a Policy1b Policy1c

%% (ii) utility-units, EZriskaversion=0 collapses Epstein-Zin to vNM
Pii=Params;
Pii.ezrisk=0;
% positive-valued utility fn
vfEZ2=struct();
vfEZ2.exoticpreferences='EpsteinZin';
vfEZ2.EZutils=1;
vfEZ2.EZpositiveutility=1;
vfEZ2.EZriskaversion='ezrisk';
vfEZ2.tolerance=vftol;
[V2a,Policy2a]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_posU,Pii,DiscountFactorParamNames,[],vfEZ2);
[V2b,Policy2b]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_posU,Pii,DiscountFactorParamNames,[],vfN);
fprintf('(ii) nod, positive utils: EZriskaversion=0 vs vNM, V, this should be zero to convergence tolerance: %2.8e \n',max(abs(V2a(:)-V2b(:))))
fprintf('(ii) nod, positive utils: EZriskaversion=0 vs vNM, raw Policy index difference, DIAGNOSTIC: %.3e at %i of %i states \n',max(abs(Policy2a(:)-Policy2b(:))),gather(sum(Policy2a(:)~=Policy2b(:))),N_a*N_z)
V2alt=ValueFnFromPolicy_InfHorz(Policy2b,0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_posU,Pii,DiscountFactorParamNames,vfEZ2);
fprintf('(ii) nod, positive utils: EZ value of the vNM policy vs V_EZ, this should be zero to convergence tolerance: %2.8e \n',max(abs(V2a(:)-V2alt(:))))
% negative-valued utility fn
vfEZ3=struct();
vfEZ3.exoticpreferences='EpsteinZin';
vfEZ3.EZutils=1;
vfEZ3.EZpositiveutility=0;
vfEZ3.EZriskaversion='ezrisk';
vfEZ3.tolerance=vftol;
[V3a,Policy3a]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_negU,Pii,DiscountFactorParamNames,[],vfEZ3);
[V3b,Policy3b]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_negU,Pii,DiscountFactorParamNames,[],vfN);
fprintf('(ii) nod, negative utils: EZriskaversion=0 vs vNM, V, this should be zero to convergence tolerance: %2.8e \n',max(abs(V3a(:)-V3b(:))))
fprintf('(ii) nod, negative utils: EZriskaversion=0 vs vNM, raw Policy index difference, DIAGNOSTIC: %.3e at %i of %i states \n',max(abs(Policy3a(:)-Policy3b(:))),gather(sum(Policy3a(:)~=Policy3b(:))),N_a*N_z)
V3alt=ValueFnFromPolicy_InfHorz(Policy3b,0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_negU,Pii,DiscountFactorParamNames,vfEZ3);
fprintf('(ii) nod, negative utils: EZ value of the vNM policy vs V_EZ, this should be zero to convergence tolerance: %2.8e \n',max(abs(V3a(:)-V3alt(:))))
% anti-vacuity: at the bank's ezrisk the EZ and vNM solves must differ
[V3d,~]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_negU,Params,DiscountFactorParamNames,[],vfEZ3);
fprintf('(ii) nod, negative utils: at ezrisk=%g vs vNM, V, should be WELL ABOVE zero: %2.8f \n',Params.ezrisk,max(abs(V3d(:)-V3b(:))))
clear V2a V2b V2alt V3a V3b V3alt V3d Policy2a Policy2b Policy3a Policy3b

%% (iv) vfoptions.EZoneminusbeta=1 versus manually scaling the return fn
Piv=Params;
Piv.ezscalefactor=(1-Params.beta)^(1/(1-1/Params.ezphi));
Piv.ezscalefactoru=1-Params.beta;
% consumption-units
ReturnFn_cons_scaled=@(aprime,a,z,r,w,ezscalefactor) ezscalefactor*EZReturnFn_cons_nod_z(aprime,a,z,r,w);
vf4=struct();
vf4.exoticpreferences='EpsteinZin';
vf4.EZutils=0;
vf4.EZriskaversion='ezgamma';
vf4.EZeis='ezphi';
vf4.tolerance=vftol;
vf4.EZoneminusbeta=1;
[V4a,Policy4a]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Piv,DiscountFactorParamNames,[],vf4);
vf4b=vf4;
vf4b.EZoneminusbeta=0;
[V4b,Policy4b]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons_scaled,Piv,DiscountFactorParamNames,[],vf4b);
fprintf('(iv) cons-units: EZoneminusbeta=1 vs manual scaling, V, this should be zero to convergence tolerance: %2.8e \n',max(abs(V4a(:)-V4b(:))))
fprintf('(iv) cons-units: EZoneminusbeta=1 vs manual scaling, raw Policy index difference, DIAGNOSTIC: %.3e at %i of %i states \n',max(abs(Policy4a(:)-Policy4b(:))),gather(sum(Policy4a(:)~=Policy4b(:))),N_a*N_z)
V4alt=ValueFnFromPolicy_InfHorz(Policy4b,0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Piv,DiscountFactorParamNames,vf4);
fprintf('(iv) cons-units: EZoneminusbeta=1 value of the scaled solve''s policy vs V, this should be zero to convergence tolerance: %2.8e \n',max(abs(V4a(:)-V4alt(:))))
% utility-units, positive
ReturnFn_posU_scaled=@(aprime,a,z,r,w,ezsigma,ezscalefactoru) ezscalefactoru*EZReturnFn_positiveUtils_nod_z(aprime,a,z,r,w,ezsigma);
vf5=struct();
vf5.exoticpreferences='EpsteinZin';
vf5.EZutils=1;
vf5.EZpositiveutility=1;
vf5.EZriskaversion='ezrisk';
vf5.tolerance=vftol;
vf5.EZoneminusbeta=1;
[V5a,Policy5a]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_posU,Piv,DiscountFactorParamNames,[],vf5);
vf5b=vf5;
vf5b.EZoneminusbeta=0;
[V5b,Policy5b]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_posU_scaled,Piv,DiscountFactorParamNames,[],vf5b);
fprintf('(iv) positive utils: EZoneminusbeta=1 vs manual scaling, V, this should be zero to convergence tolerance: %2.8e \n',max(abs(V5a(:)-V5b(:))))
fprintf('(iv) positive utils: EZoneminusbeta=1 vs manual scaling, raw Policy index difference, DIAGNOSTIC: %.3e at %i of %i states \n',max(abs(Policy5a(:)-Policy5b(:))),gather(sum(Policy5a(:)~=Policy5b(:))),N_a*N_z)
V5alt=ValueFnFromPolicy_InfHorz(Policy5b,0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_posU,Piv,DiscountFactorParamNames,vf5);
fprintf('(iv) positive utils: EZoneminusbeta=1 value of the scaled solve''s policy vs V, this should be zero to convergence tolerance: %2.8e \n',max(abs(V5a(:)-V5alt(:))))
% utility-units, negative
ReturnFn_negU_scaled=@(aprime,a,z,r,w,ezsigma,ezscalefactoru) ezscalefactoru*EZReturnFn_negativeUtils_nod_z(aprime,a,z,r,w,ezsigma);
vf6=struct();
vf6.exoticpreferences='EpsteinZin';
vf6.EZutils=1;
vf6.EZpositiveutility=0;
vf6.EZriskaversion='ezrisk';
vf6.tolerance=vftol;
vf6.EZoneminusbeta=1;
[V6a,Policy6a]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_negU,Piv,DiscountFactorParamNames,[],vf6);
vf6b=vf6;
vf6b.EZoneminusbeta=0;
[V6b,Policy6b]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_negU_scaled,Piv,DiscountFactorParamNames,[],vf6b);
fprintf('(iv) negative utils: EZoneminusbeta=1 vs manual scaling, V, this should be zero to convergence tolerance: %2.8e \n',max(abs(V6a(:)-V6b(:))))
fprintf('(iv) negative utils: EZoneminusbeta=1 vs manual scaling, raw Policy index difference, DIAGNOSTIC: %.3e at %i of %i states \n',max(abs(Policy6a(:)-Policy6b(:))),gather(sum(Policy6a(:)~=Policy6b(:))),N_a*N_z)
V6alt=ValueFnFromPolicy_InfHorz(Policy6b,0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_negU,Piv,DiscountFactorParamNames,vf6);
fprintf('(iv) negative utils: EZoneminusbeta=1 value of the scaled solve''s policy vs V, this should be zero to convergence tolerance: %2.8e \n',max(abs(V6a(:)-V6alt(:))))
% anti-vacuity: EZoneminusbeta=1 must change the answer relative to EZoneminusbeta=0 with the UNscaled return
[V4c,~]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Piv,DiscountFactorParamNames,[],vf4b);
fprintf('(iv) cons-units: EZoneminusbeta=1 vs =0 on the same return fn, V, should be WELL ABOVE zero: %2.8f \n',max(abs(V4a(:)-V4c(:))))
clear V4a V4b V4c V4alt V5a V5b V5alt V6a V6b V6alt Policy4a Policy4b Policy5a Policy5b Policy6a Policy6b

%% (xii) Homogeneity (units invariance) in cons-units, plus a high-risk-aversion stress leg
% Traditional Epstein-Zin (EZutils=0) is homogeneous of degree 1 in consumption: scaling the
% units of the model by lambda (a_grid*lambda and the level parameter w; r and z are unit-free)
% scales c by lambda everywhere, so V must equal lambda*V_baseline (checked in relative terms)
% with Policy indices identical (the identical computation up to scale, so ties scale too).
% Utility-units EZ is NOT homogeneous of degree 1, so this is cons-units only.
vf12=struct();
vf12.exoticpreferences='EpsteinZin';
vf12.EZutils=0;
vf12.EZriskaversion='ezgamma';
vf12.EZeis='ezphi';
vf12.tolerance=vftol;
for gg=1:2
    P12=Params;
    if gg==1
        gammastr=sprintf('ezgamma=%g',Params.ezgamma); % baseline risk aversion
    else
        P12.ezgamma=20; % high-risk-aversion overflow stress
        gammastr='ezgamma=20';
    end
    [Vbase,Policybase]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,P12,DiscountFactorParamNames,[],vf12);
    fprintf('(xii) homogeneity leg (%s), all(isfinite(V)), this should be one: %i \n',gammastr,gather(all(isfinite(Vbase(:)))))
    for lam=[10^(-3),10^3]
        Plam=P12;
        Plam.w=lam*P12.w;
        % The solver's tolerance is ABSOLUTE on V, and V scales by lambda, so the scaled solve gets
        % tolerance lambda*vftol: the same RELATIVE accuracy as the baseline, which is what the
        % relative check below asks for. At lambda=1e-3 and the unscaled 1e-10 the fixed point was
        % only pinned to tolerance/(1-beta)/|V| = 4e-7 relative, and the ezgamma=20 leg read
        % 5.19e-7 against the 1e-7 bar on 2026-10-27 (the ezgamma=3 leg happened to land well
        % inside its tolerance). Units-consistent tolerance, not a looser bar.
        vf12lam=vf12;
        vf12lam.tolerance=lam*vftol;
        [Vlam,Policylam]=ValueFnIter_InfHorz(0,n_a,n_z,[],lam*a_grid,z_grid,pi_z,ReturnFn_cons,Plam,DiscountFactorParamNames,[],vf12lam);
        fprintf('(xii) homogeneity, lambda=%g (%s), V relative, this should be zero to convergence tolerance: %2.8e \n',lam,gammastr,max(abs(Vlam(:)-lam*Vbase(:)))/max(abs(Vlam(:))))
        fprintf('(xii) homogeneity, lambda=%g (%s), Policy, this should be zero: %.3e \n',lam,gammastr,max(abs(Policylam(:)-Policybase(:))))
        fprintf('(xii) homogeneity, lambda=%g (%s), all(isfinite(V)), this should be one: %i \n',lam,gammastr,gather(all(isfinite(Vlam(:)))))
    end
end
clear Vbase Policybase Vlam Policylam

%% (i) and (ii) again with d (ValueFnIter_InfHorz_EpsteinZin_raw is a different file from the nod raw)
fprintf('--- with d --- \n')
[V1a,Policy1a]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons_d,Pi,DiscountFactorParamNames,[],vfEZ1);
[V1b,Policy1b]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU_d,Pi,DiscountFactorParamNames,[],vfN);
V1bt=((1-Pi.ezgamma)*V1b).^(1/(1-Pi.ezgamma));
fprintf('(i) d: gamma=1/phi collapse, V after transform (relative), this should be zero to convergence tolerance: %2.8e \n',max(abs(V1a(:)-V1bt(:)))/max(abs(V1a(:))))
fprintf('(i) d: gamma=1/phi collapse, raw (d,aprime) index difference, DIAGNOSTIC: %.3e at %i of %i state-channels \n',max(abs(Policy1a(:)-Policy1b(:))),gather(sum(Policy1a(:)~=Policy1b(:))),2*N_a*N_z)
V1alt=ValueFnFromPolicy_InfHorz(Policy1b,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons_d,Pi,DiscountFactorParamNames,vfEZ1);
fprintf('(i) d: gamma=1/phi collapse, EZ value of the vNM policy vs V_EZ (relative), this should be zero to convergence tolerance: %2.8e \n',max(abs(V1a(:)-V1alt(:)))/max(abs(V1a(:))))
[V2a,Policy2a]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU_d,Pii,DiscountFactorParamNames,[],vfEZ2);
[V2b,Policy2b]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU_d,Pii,DiscountFactorParamNames,[],vfN);
fprintf('(ii) d, positive utils: EZriskaversion=0 vs vNM, V, this should be zero to convergence tolerance: %2.8e \n',max(abs(V2a(:)-V2b(:))))
fprintf('(ii) d, positive utils: EZriskaversion=0 vs vNM, raw (d,aprime) index difference, DIAGNOSTIC: %.3e at %i of %i state-channels \n',max(abs(Policy2a(:)-Policy2b(:))),gather(sum(Policy2a(:)~=Policy2b(:))),2*N_a*N_z)
V2alt=ValueFnFromPolicy_InfHorz(Policy2b,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_posU_d,Pii,DiscountFactorParamNames,vfEZ2);
fprintf('(ii) d, positive utils: EZ value of the vNM policy vs V_EZ, this should be zero to convergence tolerance: %2.8e \n',max(abs(V2a(:)-V2alt(:))))
[V3a,Policy3a]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU_d,Pii,DiscountFactorParamNames,[],vfEZ3);
[V3b,Policy3b]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU_d,Pii,DiscountFactorParamNames,[],vfN);
fprintf('(ii) d, negative utils: EZriskaversion=0 vs vNM, V, this should be zero to convergence tolerance: %2.8e \n',max(abs(V3a(:)-V3b(:))))
fprintf('(ii) d, negative utils: EZriskaversion=0 vs vNM, raw (d,aprime) index difference, DIAGNOSTIC: %.3e at %i of %i state-channels \n',max(abs(Policy3a(:)-Policy3b(:))),gather(sum(Policy3a(:)~=Policy3b(:))),2*N_a*N_z)
V3alt=ValueFnFromPolicy_InfHorz(Policy3b,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU_d,Pii,DiscountFactorParamNames,vfEZ3);
fprintf('(ii) d, negative utils: EZ value of the vNM policy vs V_EZ, this should be zero to convergence tolerance: %2.8e \n',max(abs(V3a(:)-V3alt(:))))

%%
output=struct();

end
