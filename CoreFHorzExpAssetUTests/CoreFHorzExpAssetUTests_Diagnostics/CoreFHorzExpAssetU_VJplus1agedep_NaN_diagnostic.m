% Diagnostic for the 12 'V_Jplus1 with age-dependent shocks (jstar=7), this should be zero: NaN'
% checks in CoreFHorzExpAssetUTests parts 1 and 3 (the noa1 tiers), seen on the 2026-10-12 run.
%
% Why this needs a run: MATLAB's max() OMITS NaN, so the printed NaN means EVERY element of
% abs(Vbase(:)-Vshort(:)) is NaN, not merely the infeasible a=0 states. That rules out the easy
% explanation (-Inf minus -Inf at a=0 only) and points at one of the two arrays being NaN
% throughout. This script says which, and where.
%
% Mirrors CoreFHorzExpAssetU_nod1_z_noe_noa1_nosemiz, which is one of the 6 failing variants.
% Run from the CoreFHorzExpAssetUTests folder.

% Re-running the setup is harmless (it only redefines grids/params) and keeps this script
% runnable on its own; when called from CoreFHorzExpAssetUTests the bank has already run it.
addpath('./CoreFHorzExpAssetUTests_Setup/')
addpath('./CoreFHorzExpAssetU_ReturnFns/Noa1_ReturnFns/')
CoreFHorzExpAssetU_setup

n_d=n_d_withoutd1; d_grid=d_grid_withoutd1;
n_a=n_a_justexpasset; a_grid=a_grid_justexpasset;

vfoptions=struct();
vfoptions.experienceassetu=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions.n_u=vfoptionsbaseline.n_u;
vfoptions.u_grid=vfoptionsbaseline.u_grid;
vfoptions.pi_u=vfoptionsbaseline.pi_u;

ReturnFn=@(d2,a,z,r,w,kappa_j,sigma,agej,Jr,pension,transfer) ReturnFn_nod1_z_noe_noa1_nosemiz(d2,a,z,r,w,kappa_j,sigma,agej,Jr,pension,transfer);

%% Reproduce the age-dependent-shocks V_Jplus1 test exactly as the subcode does
jstar=round(N_j/3);
Njs=jstar-1;
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
pi_z_J=pi_z.*ones(1,1,N_j);
pi_z_J(:,:,1:2:N_j)=0.5*pi_z_J(:,:,1:2:N_j)+0.5*eye(n_z);
vfoptionsjs=vfoptions;
[Vbase,Policybase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
vfoptionsjs.V_Jplus1=Vbase(:,:,jstar);
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
[Vshort,Policyshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z_J(:,:,1:Njs),ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);

fprintf('\n================ V_Jplus1 age-dependent NaN diagnostic ================\n')
fprintf('N_j=%i  jstar=%i  Njs=%i   n_a=%i  n_z=%i  n_u=%i\n',N_j,jstar,Njs,n_a,n_z,vfoptions.n_u)
fprintf('Jr=%i   kappa_j(1:%i)=[%s]\n',Params.Jr,Njs,num2str(Params.kappa_j(1:Njs),'%.2f '))
fprintf('a_grid(1)=%g (0 makes working-age c=w*kappa_j*d2*a*z exactly 0, so F=-Inf there)\n',a_grid(1))

fprintf('\n-- the printed check --\n')
fprintf('  max(abs(Vbase-Vshort))      = %.3e   (NaN here means EVERY element is NaN)\n',max(abs(Vbase(:)-Vshort(:))))
fprintf('  max(abs(Policybase-Policyshort)) = %.3e\n',max(abs(Policybase(:)-Policyshort(:))))

fprintf('\n-- where the badness is --\n')
for nm={'Vbase','Vshort'}
    X=eval(nm{1});
    fprintf('  %-7s numel=%7i   NaN=%7i   -Inf=%7i   +Inf=%7i   finite=%7i\n', ...
        nm{1},numel(X),sum(isnan(X(:))),sum(X(:)==-Inf),sum(X(:)==Inf),sum(isfinite(X(:))));
end
D=abs(Vbase(:)-Vshort(:));
fprintf('  diff:   numel=%7i   NaN=%7i   finite=%7i\n',numel(D),sum(isnan(D)),sum(isfinite(D)));
if any(isfinite(D))
    fprintf('  max over the FINITE diffs only = %.3e   <-- this is the number the check meant to print\n',max(D(isfinite(D))));
end

fprintf('\n-- NaN by age (j), to see whether it is all ages or just some --\n')
for jj=1:Njs
    vb=Vbase(:,:,jj); vs=Vshort(:,:,jj);
    fprintf('  j=%2i  Vbase: NaN=%4i -Inf=%4i finite=%4i |  Vshort: NaN=%4i -Inf=%4i finite=%4i\n', ...
        jj,sum(isnan(vb(:))),sum(vb(:)==-Inf),sum(isfinite(vb(:))), ...
        sum(isnan(vs(:))),sum(vs(:)==-Inf),sum(isfinite(vs(:))));
end

fprintf('\n-- a=1 (the a=0 grid point) vs a>1, at j=1 --\n')
fprintf('  Vbase(a=1,:,1)  = [%s]\n',num2str(Vbase(1,:,1),'%.4g '))
fprintf('  Vshort(a=1,:,1) = [%s]\n',num2str(Vshort(1,:,1),'%.4g '))
fprintf('  Vbase(a=2,:,1)  = [%s]\n',num2str(Vbase(2,:,1),'%.4g '))
fprintf('  Vshort(a=2,:,1) = [%s]\n',num2str(Vshort(2,:,1),'%.4g '))
fprintf('\nVERDICT GUIDE: if both arrays are -Inf at the SAME places and finite elsewhere with the\n')
fprintf('finite diffs at 0, the check is a comparison artefact (-Inf minus -Inf) and the bank line\n')
fprintf('should mask non-finite entries. If instead one array is NaN where the other is finite,\n')
fprintf('that is a real defect in the V_Jplus1 branch of the ExpAssetu noa1 raws.\n')
fprintf('=======================================================================\n')

fprintf('DIAGNOSTIC COMPLETE\n')
