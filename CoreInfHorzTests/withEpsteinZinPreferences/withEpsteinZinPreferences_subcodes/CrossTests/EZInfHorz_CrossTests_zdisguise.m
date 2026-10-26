function output=EZInfHorz_CrossTests_zdisguise(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames)
% Cross-tests: run the same model disguised as a different z process, and check the answer follows.
%
% A: SPLIT. Replace each z state by two identical copies, and split the transition probabilities
%    evenly between them. That is the very same stochastic process, so V and Policy must simply
%    be duplicated column by column. This is the check that would catch a raw that indexed pi_z
%    by the wrong dimension, or summed over z rather than z', inside the certainty equivalent.
%
% B: RELABEL. Permute the z states (grid and transition matrix together). V and Policy must
%    permute with them. A raw that mixed up z and zprime anywhere passes A but fails B, because
%    a permuted pi_z is no longer symmetric in the way the split one is.
%
% C: DEGENERATE. n_z points all equal to z=1 with uniform pi_z: the solution must be identical
%    across the z slices AND reproduce a 1-point-z solve (n_z=1, z_grid=1, pi_z=1). Each
%    degenerate solve still exercises the full certainty-equivalent machinery (the expectation is
%    over five equal values). This is also the nearest thing to a no-shock model the bank can
%    run, since n_z=0 is refused (see the guards cross-test).
%
% D: TWO-DIMENSIONAL z. Two independent z variables entering the return fn only through their
%    product z1*z2 is the SAME model as one z whose grid holds the products in the joint order,
%    with the joint transition matrix. The z_gridvals defect that the QH bank found in 87bca184
%    (the exotic-preference branches were handed z_grid where every other branch gets z_gridvals,
%    so a two-dimensional z read column 2 of a one-column array) had the same shape in the
%    Epstein-Zin branch and was fixed with it, with no Epstein-Zin coverage until now. Run on
%    both lowmemory paths: the lowmemory raws convert z to gridvals themselves, a separate path.
%    With a discriminator: z2 spread vs z2 degenerate must differ, or the test could pass by z2
%    being silently dropped.
%
% A, B, C and D run for all three EZ cases without d; A and B also with d for cons-units. All
% are the identical computation up to floating-point order (same raw, same operator), so Policy
% is graded on the raw indices, as the QH bank does for these.

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn_cons=@(aprime,a,z,r,w) EZReturnFn_cons_nod_z(aprime,a,z,r,w);
ReturnFn_posU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_positiveUtils_nod_z(aprime,a,z,r,w,ezsigma);
ReturnFn_negU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_negativeUtils_nod_z(aprime,a,z,r,w,ezsigma);
ReturnFn_cons_d=@(d,aprime,a,z,r,w,varphi) EZReturnFn_cons_d_z(d,aprime,a,z,r,w,varphi);
% two-dimensional z versions (z1*z2 in place of z)
ReturnFn_cons_2z=@(aprime,a,z1,z2,r,w) EZReturnFn_cons_nod_z(aprime,a,z1*z2,r,w);
ReturnFn_posU_2z=@(aprime,a,z1,z2,r,w,ezsigma) EZReturnFn_positiveUtils_nod_z(aprime,a,z1*z2,r,w,ezsigma);
ReturnFn_negU_2z=@(aprime,a,z1,z2,r,w,ezsigma) EZReturnFn_negativeUtils_nod_z(aprime,a,z1*z2,r,w,ezsigma);

vftol=10^(-10); % see the driver header

fprintf('\n===== z-disguise cross-tests ===== \n')

% A: each z state split into two identical copies, transitions split evenly between the copies
n_z_split=2*N_z;
z_grid_split=repelem(z_grid,2,1);
pi_z_split=repelem(pi_z,2,2)/2;
fprintf('split pi_z rows sum to one, this should be zero: %2.8e \n',max(abs(sum(pi_z_split,2)-1)))

% B: relabelling, interleave-split so the permuted grid is clearly not monotone
perm=[(2:2:N_z)';(1:2:N_z)'];
z_grid_perm=z_grid(perm);
pi_z_perm=pi_z(perm,perm);
fprintf('permuted pi_z rows sum to one, this should be zero: %2.8e \n',max(abs(sum(pi_z_perm,2)-1)))

% C: degenerate z
zdeg_grid=ones(N_z,1);
pi_zdeg=ones(N_z,N_z)/N_z;

% D: two-dimensional z. z1 is the FASTEST-moving index in the joint state, so pi_z2 is the left
% kron factor. Get the ordering wrong and the two models below are genuinely different.
n_z_2d=[3,2];
z1_grid=z_grid(1:3);
z2_grid=[0.9;1.1];
z_grid_2d=[z1_grid;z2_grid];
pi_z1=pi_z(1:3,1:3)./sum(pi_z(1:3,1:3),2);
pi_z2=[0.7,0.3;0.3,0.7];
pi_z_2d=kron(pi_z2,pi_z1);
zprod_grid=kron(z2_grid,ones(3,1)).*kron(ones(2,1),z1_grid); % the products, in the same joint order (z1 fastest)

for ezcase=1:3
    vfoptions=struct();
    vfoptions.exoticpreferences='EpsteinZin';
    vfoptions.tolerance=vftol;
    if ezcase==1
        casestr='cons-units';
        vfoptions.EZutils=0;
        vfoptions.EZriskaversion='ezgamma';
        vfoptions.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
        ReturnFn_2z=ReturnFn_cons_2z;
    elseif ezcase==2
        casestr='positive utils';
        vfoptions.EZutils=1;
        vfoptions.EZpositiveutility=1;
        vfoptions.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
        ReturnFn_2z=ReturnFn_posU_2z;
    else
        casestr='negative utils';
        vfoptions.EZutils=1;
        vfoptions.EZpositiveutility=0;
        vfoptions.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
        ReturnFn_2z=ReturnFn_negU_2z;
    end
    fprintf('--- nod [EZ %s] --- \n',casestr)

    %% A and B
    [V0,P0]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    [Vsp,Psp]=ValueFnIter_InfHorz(0,n_a,n_z_split,[],a_grid,z_grid_split,pi_z_split,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    fprintf('A nod [EZ %s]: split V vs duplicated base V, this should be zero: %2.8e \n',casestr,max(abs(Vsp(:)-reshape(repelem(V0,1,2),[N_a*n_z_split,1]))))
    fprintf('A nod [EZ %s]: split Policy vs duplicated base Policy, this should be zero: %.3e \n',casestr,max(abs(Psp(:)-reshape(repelem(P0,1,1,2),[N_a*n_z_split,1]))))

    [Vpm,Ppm]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid_perm,pi_z_perm,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    fprintf('B nod [EZ %s]: permuted V vs permuted base V, this should be zero: %2.8e \n',casestr,max(abs(Vpm(:)-reshape(V0(:,perm),[N_a*N_z,1]))))
    fprintf('B nod [EZ %s]: permuted Policy vs permuted base Policy, this should be zero: %.3e \n',casestr,max(abs(Ppm(:)-reshape(P0(:,:,perm),[N_a*N_z,1]))))

    %% C: degenerate z vs a 1-point z
    [Vdeg,Pdeg]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,zdeg_grid,pi_zdeg,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    Vslicedev=Vdeg-Vdeg(:,1); % implicit expansion: every z slice against the first
    Pslicedev=Pdeg-Pdeg(:,:,1);
    fprintf('C nod [EZ %s]: degenerate z, V identical across z slices, this should be zero: %2.8e \n',casestr,max(abs(Vslicedev(:))))
    fprintf('C nod [EZ %s]: degenerate z, Policy identical across z slices, this should be zero: %.3e \n',casestr,max(abs(Pslicedev(:))))
    [Vonept,Ponept]=ValueFnIter_InfHorz(0,n_a,1,[],a_grid,1,1,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    fprintf('C nod [EZ %s]: degenerate z vs 1-point z, V, this should be zero: %2.8e \n',casestr,max(abs(Vdeg(:,1)-Vonept(:))))
    fprintf('C nod [EZ %s]: degenerate z vs 1-point z, Policy, this should be zero: %.3e \n',casestr,max(abs(reshape(Pdeg(:,:,1),[N_a,1])-reshape(Ponept,[N_a,1]))))
    % and the degenerate model must differ from the risky one (else C is a no-shock model twice over)
    fprintf('C nod [EZ %s]: degenerate z vs the bank''s z, V, should be WELL ABOVE zero: %2.8f \n',casestr,max(abs(Vdeg(:)-V0(:))))

    %% D: two-dimensional z vs the equivalent single z, on both lowmemory paths
    for lm=0:1
        vf2d=vfoptions;
        vf2d.lowmemory=lm;
        [V2d,P2d]=ValueFnIter_InfHorz(0,n_a,n_z_2d,[],a_grid,z_grid_2d,pi_z_2d,ReturnFn_2z,Params,DiscountFactorParamNames,[],vf2d);
        [V2done,P2done]=ValueFnIter_InfHorz(0,n_a,prod(n_z_2d),[],a_grid,zprod_grid,pi_z_2d,ReturnFn,Params,DiscountFactorParamNames,[],vf2d);
        fprintf('D nod, lowmemory=%i [EZ %s]: two-dimensional n_z vs the equivalent single z, V, this should be zero: %2.8e \n',lm,casestr,max(abs(reshape(V2d,[N_a*prod(n_z_2d),1])-reshape(V2done,[N_a*prod(n_z_2d),1]))))
        fprintf('D nod, lowmemory=%i [EZ %s]: two-dimensional n_z vs the equivalent single z, Policy, this should be zero: %.3e \n',lm,casestr,max(abs(double(reshape(P2d,[],1))-double(reshape(P2done,[],1)))))
        % Discriminator: the two-dimensional solve must actually depend on z2
        [V2dflat,~]=ValueFnIter_InfHorz(0,n_a,n_z_2d,[],a_grid,[z1_grid;[1;1]],pi_z_2d,ReturnFn_2z,Params,DiscountFactorParamNames,[],vf2d);
        fprintf('D nod, lowmemory=%i [EZ %s]: z2 actually reaches the return fn, max|V(z2 spread) - V(z2 degenerate)| should NOT be zero: %2.8f \n',lm,casestr,max(abs(reshape(V2d,[N_a*prod(n_z_2d),1])-reshape(V2dflat,[N_a*prod(n_z_2d),1]))))
    end

    clear V0 P0 Vsp Psp Vpm Ppm Vdeg Pdeg Vonept Ponept V2d P2d V2done P2done V2dflat
end

%% A and B with d, cons-units
fprintf('--- d [EZ cons-units] --- \n')
vfoptions=struct();
vfoptions.exoticpreferences='EpsteinZin';
vfoptions.tolerance=vftol;
vfoptions.EZutils=0;
vfoptions.EZriskaversion='ezgamma';
vfoptions.EZeis='ezphi';
[V0,P0]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons_d,Params,DiscountFactorParamNames,[],vfoptions);
[Vsp,Psp]=ValueFnIter_InfHorz(n_d,n_a,n_z_split,d_grid,a_grid,z_grid_split,pi_z_split,ReturnFn_cons_d,Params,DiscountFactorParamNames,[],vfoptions);
fprintf('A d [EZ cons-units]: split V vs duplicated base V, this should be zero: %2.8e \n',max(abs(Vsp(:)-reshape(repelem(V0,1,2),[N_a*n_z_split,1]))))
fprintf('A d [EZ cons-units]: split Policy vs duplicated base Policy, this should be zero: %.3e \n',max(abs(Psp(:)-reshape(repelem(P0,1,1,2),[2*N_a*n_z_split,1]))))
[Vpm,Ppm]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid_perm,pi_z_perm,ReturnFn_cons_d,Params,DiscountFactorParamNames,[],vfoptions);
fprintf('B d [EZ cons-units]: permuted V vs permuted base V, this should be zero: %2.8e \n',max(abs(Vpm(:)-reshape(V0(:,perm),[N_a*N_z,1]))))
fprintf('B d [EZ cons-units]: permuted Policy vs permuted base Policy, this should be zero: %.3e \n',max(abs(Ppm(:)-reshape(P0(:,:,perm),[2*N_a*N_z,1]))))

%%
output=struct();

end
