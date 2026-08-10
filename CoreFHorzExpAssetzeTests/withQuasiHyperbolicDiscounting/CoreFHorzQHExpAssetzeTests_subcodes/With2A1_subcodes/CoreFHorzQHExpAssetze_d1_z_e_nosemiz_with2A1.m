function output=CoreFHorzQHExpAssetze_d1_z_e_nosemiz_with2A1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% Quasi-hyperbolic experienceassetze with TWO standard endogenous assets -> triggers the
% DC2A / GI2A / DC2A_GI2A QH-ExpAssetze code paths (length(n_a1)>1). a1_1 is divide-conquered,
% the binary a1_2 is folded, a2 is the experience asset. Runs Naive then Sophisticated, with a
% continuation-value (Valt) check beside every V check, then the exponential cross-tests.
% shocks: {z (markov), e (iid)} -> valid lowmemory {0,1,2}.

% Build the binary second standard endogenous asset a1_2, inserted between a1_1 and a2
n_a1_1=n_a(1); n_a2exp=n_a(2);
a1_1_grid=a_grid(1:n_a1_1);
a2_grid=a_grid(n_a1_1+1:end);
a1_2_grid=[0;1]; % binary second asset (capped high-return asset)
n_a=[n_a1_1,2,n_a2exp];
a_grid=[a1_1_grid;a1_2_grid;a2_grid];
Params.r2=0.08; % return on the binary asset (higher than r, so it is used up to the cap)

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.pi_e=simoptionsbaseline.pi_e;
simoptions.e_grid=simoptionsbaseline.e_grid;

ReturnFn=@(d1,d2,a1prime,a1_2prime,a1,a1_2,a2,z,e,r,r2,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_ExpAssetze_d1_z_e_with2A1(d1,d2,a1prime,a1_2prime,a1,a1_2,a2,z,e,r,r2,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);

% Experience asset (z variant)
vfoptions.experienceassetze=1;
simoptions.experienceassetze=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;
simoptions.z_grid=z_grid;

%% Quasi-Hyperbolic Discounting
vfoptions.exoticpreferences='QuasiHyperbolic';
vfoptions.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;

for qhcase=1:2
    if qhcase==1
        vfoptions.quasi_hyperbolic='Naive'; qh='Naive';
    else
        vfoptions.quasi_hyperbolic='Sophisticated'; qh='Sophisticated';
    end

    % Base (2 standard assets -> combined N_a1)
    vfoptions1=vfoptions;
    if qhcase==1 % Naive also returns Policyalt (4th output), needed by the ValueFnFromPolicy oracle below
        [V1,Policy1,V1alt,Policy1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    else
        [V1,Policy1,V1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    end

    % ValueFnFromPolicy oracle on the base method (exercises QH+experienceassetze reconstruction with two standard assets; GI2A/DC2A_GI2A not checked -- QH ValueFnFromPolicy has no gridinterplayer variant yet)
    vfoptionsVFP=vfoptions1; vfoptionsVFP.lowmemory=0;
    if qhcase==1, vfoptionsVFP.Policyalt=Policy1alt; end
    [V1fromPolicy,V1altfromPolicy]=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP);
    fprintf('%s ValueFnFromPolicy, this should be zero: %2.8f \n',qh,max(abs(V1fromPolicy(:)-V1(:))))
    fprintf('%s ValueFnFromPolicy (Valt), this should be zero: %2.8f \n',qh,max(abs(V1altfromPolicy(:)-V1alt(:))))

    % Divide-and-conquer -> DC2A, should give same answer
    vfoptions2=vfoptions; vfoptions2.divideandconquer=1;
    if qhcase==1 % Naive also returns Policyalt (4th output), needed by the Policyalt equivalence checks below
        [V2,Policy2,V2alt,Policy2alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    else
        [V2,Policy2,V2alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    end
    fprintf('%s DC2A, this should be zero: %2.8f \n',qh,max(abs(V1(:)-V2(:))))
    fprintf('%s DC2A (Valt), this should be zero: %2.8f \n',qh,max(abs(V1alt(:)-V2alt(:))))
    fprintf('%s DC2A (Policy), this should be zero: %2.8f \n',qh,max(abs(Policy1(:)-Policy2(:))))
    if qhcase==1, fprintf('%s DC2A (Policyalt), this should be zero: %2.8f \n',qh,max(abs(Policy1alt(:)-Policy2alt(:)))); end

    % lowmemory on base
    vfoptions1.lowmemory=1;
    if qhcase==1
        [V1B,Policy1B,V1Balt,Policy1Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    else
        [V1B,Policy1B,V1Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    end
    fprintf('%s lowmemory=1, this should be zero: %2.8f \n',qh,max(abs(V1(:)-V1B(:))))
    fprintf('%s lowmemory=1 (Valt), this should be zero: %2.8f \n',qh,max(abs(V1alt(:)-V1Balt(:))))
    fprintf('%s lowmemory=1 (Policy), this should be zero: %2.8f \n',qh,max(abs(Policy1(:)-Policy1B(:))))
    if qhcase==1, fprintf('%s lowmemory=1 (Policyalt), this should be zero: %2.8f \n',qh,max(abs(Policy1alt(:)-Policy1Balt(:)))); end
    vfoptions1.lowmemory=2;
    if qhcase==1
        [V1C,Policy1C,V1Calt,Policy1Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    else
        [V1C,Policy1C,V1Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
    end
    fprintf('%s lowmemory=2, this should be zero: %2.8f \n',qh,max(abs(V1(:)-V1C(:))))
    fprintf('%s lowmemory=2 (Valt), this should be zero: %2.8f \n',qh,max(abs(V1alt(:)-V1Calt(:))))
    fprintf('%s lowmemory=2 (Policy), this should be zero: %2.8f \n',qh,max(abs(Policy1(:)-Policy1C(:))))
    if qhcase==1, fprintf('%s lowmemory=2 (Policyalt), this should be zero: %2.8f \n',qh,max(abs(Policy1alt(:)-Policy1Calt(:)))); end
    vfoptions1.lowmemory=0;
    % lowmemory on DC2A
    vfoptions2.lowmemory=1;
    if qhcase==1
        [V2B,Policy2B,V2Balt,Policy2Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    else
        [V2B,Policy2B,V2Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    end
    fprintf('%s lowmemory=1 (DC2A), this should be zero: %2.8f \n',qh,max(abs(V2(:)-V2B(:))))
    fprintf('%s lowmemory=1 (DC2A, Valt), this should be zero: %2.8f \n',qh,max(abs(V2alt(:)-V2Balt(:))))
    fprintf('%s lowmemory=1 (DC2A, Policy), this should be zero: %2.8f \n',qh,max(abs(Policy2(:)-Policy2B(:))))
    if qhcase==1, fprintf('%s lowmemory=1 (DC2A, Policyalt), this should be zero: %2.8f \n',qh,max(abs(Policy2alt(:)-Policy2Balt(:)))); end
    vfoptions2.lowmemory=2;
    if qhcase==1
        [V2C,Policy2C,V2Calt,Policy2Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    else
        [V2C,Policy2C,V2Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
    end
    fprintf('%s lowmemory=2 (DC2A), this should be zero: %2.8f \n',qh,max(abs(V2(:)-V2C(:))))
    fprintf('%s lowmemory=2 (DC2A, Valt), this should be zero: %2.8f \n',qh,max(abs(V2alt(:)-V2Calt(:))))
    fprintf('%s lowmemory=2 (DC2A, Policy), this should be zero: %2.8f \n',qh,max(abs(Policy2(:)-Policy2C(:))))
    if qhcase==1, fprintf('%s lowmemory=2 (DC2A, Policyalt), this should be zero: %2.8f \n',qh,max(abs(Policy2alt(:)-Policy2Calt(:)))); end
    vfoptions2.lowmemory=0;

    % Grid interpolation -> GI2A
    vfoptions3=vfoptions; vfoptions3.gridinterplayer=1; vfoptions3.ngridinterp=5;
    if qhcase==1 % Naive also returns Policyalt (4th output), needed by the ValueFnFromPolicy oracle below
        [V3,Policy3,V3alt,Policy3alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    else
        [V3,Policy3,V3alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    end
    % ValueFnFromPolicy oracle on GI2A
    vfoptionsVFP3=vfoptions3; vfoptionsVFP3.lowmemory=0;
    if qhcase==1, vfoptionsVFP3.Policyalt=Policy3alt; end
    [V3fromPolicy,V3altfromPolicy]=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptionsVFP3);
    fprintf('%s ValueFnFromPolicy (GI2A), this should be zero: %2.8f \n',qh,max(abs(V3fromPolicy(:)-V3(:))))
    fprintf('%s ValueFnFromPolicy (GI2A, Valt), this should be zero: %2.8f \n',qh,max(abs(V3altfromPolicy(:)-V3alt(:))))
    % DC + GI -> DC2A_GI2A, should match GI2A
    vfoptions4=vfoptions; vfoptions4.divideandconquer=1; vfoptions4.gridinterplayer=1; vfoptions4.ngridinterp=5;
    if qhcase==1 % Naive also returns Policyalt (4th output), needed by the Policyalt equivalence checks below
        [V4,Policy4,V4alt,Policy4alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    else
        [V4,Policy4,V4alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    end
    fprintf('%s DC2A_GI2A vs GI2A, this should be zero: %2.8f \n',qh,max(abs(V3(:)-V4(:))))
    fprintf('%s DC2A_GI2A vs GI2A (Valt), this should be zero: %2.8f \n',qh,max(abs(V3alt(:)-V4alt(:))))
    fprintf('%s DC2A_GI2A vs GI2A (Policy), this should be zero: %2.8f \n',qh,max(abs(Policy3(:)-Policy4(:))))
    if qhcase==1, fprintf('%s DC2A_GI2A vs GI2A (Policyalt), this should be zero: %2.8f \n',qh,max(abs(Policy3alt(:)-Policy4alt(:)))); end

    % lowmemory on GI2A
    vfoptions3.lowmemory=1;
    if qhcase==1
        [V3B,Policy3B,V3Balt,Policy3Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    else
        [V3B,Policy3B,V3Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    end
    fprintf('%s lowmemory=1 (GI2A), this should be zero: %2.8f \n',qh,max(abs(V3(:)-V3B(:))))
    fprintf('%s lowmemory=1 (GI2A, Valt), this should be zero: %2.8f \n',qh,max(abs(V3alt(:)-V3Balt(:))))
    fprintf('%s lowmemory=1 (GI2A, Policy), this should be zero: %2.8f \n',qh,max(abs(Policy3(:)-Policy3B(:))))
    if qhcase==1, fprintf('%s lowmemory=1 (GI2A, Policyalt), this should be zero: %2.8f \n',qh,max(abs(Policy3alt(:)-Policy3Balt(:)))); end
    vfoptions3.lowmemory=2;
    if qhcase==1
        [V3C,Policy3C,V3Calt,Policy3Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    else
        [V3C,Policy3C,V3Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
    end
    fprintf('%s lowmemory=2 (GI2A), this should be zero: %2.8f \n',qh,max(abs(V3(:)-V3C(:))))
    fprintf('%s lowmemory=2 (GI2A, Valt), this should be zero: %2.8f \n',qh,max(abs(V3alt(:)-V3Calt(:))))
    fprintf('%s lowmemory=2 (GI2A, Policy), this should be zero: %2.8f \n',qh,max(abs(Policy3(:)-Policy3C(:))))
    if qhcase==1, fprintf('%s lowmemory=2 (GI2A, Policyalt), this should be zero: %2.8f \n',qh,max(abs(Policy3alt(:)-Policy3Calt(:)))); end
    vfoptions3.lowmemory=0;
    % lowmemory on DC2A_GI2A
    vfoptions4.lowmemory=1;
    if qhcase==1
        [V4B,Policy4B,V4Balt,Policy4Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    else
        [V4B,Policy4B,V4Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    end
    fprintf('%s lowmemory=1 (DC2A_GI2A), this should be zero: %2.8f \n',qh,max(abs(V4(:)-V4B(:))))
    fprintf('%s lowmemory=1 (DC2A_GI2A, Valt), this should be zero: %2.8f \n',qh,max(abs(V4alt(:)-V4Balt(:))))
    fprintf('%s lowmemory=1 (DC2A_GI2A, Policy), this should be zero: %2.8f \n',qh,max(abs(Policy4(:)-Policy4B(:))))
    if qhcase==1, fprintf('%s lowmemory=1 (DC2A_GI2A, Policyalt), this should be zero: %2.8f \n',qh,max(abs(Policy4alt(:)-Policy4Balt(:)))); end
    vfoptions4.lowmemory=2;
    if qhcase==1
        [V4C,Policy4C,V4Calt,Policy4Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    else
        [V4C,Policy4C,V4Calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
    end
    fprintf('%s lowmemory=2 (DC2A_GI2A), this should be zero: %2.8f \n',qh,max(abs(V4(:)-V4C(:))))
    fprintf('%s lowmemory=2 (DC2A_GI2A, Valt), this should be zero: %2.8f \n',qh,max(abs(V4alt(:)-V4Calt(:))))
    fprintf('%s lowmemory=2 (DC2A_GI2A, Policy), this should be zero: %2.8f \n',qh,max(abs(Policy4(:)-Policy4C(:))))
    if qhcase==1, fprintf('%s lowmemory=2 (DC2A_GI2A, Policyalt), this should be zero: %2.8f \n',qh,max(abs(Policy4alt(:)-Policy4Calt(:)))); end
    vfoptions4.lowmemory=0;
end

%% Versus exponential discounting (at baseline; DC2A already shown equal to base above)
% (i) at the actual QH beta0, Naive's continuation value equals the exponential value function
vfoptionsE=vfoptions; vfoptionsE.exoticpreferences='None';
[Vexp,~]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsE);
vfoptionsN=vfoptions; vfoptionsN.quasi_hyperbolic='Naive';
[~,~,VnAalt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsN);
fprintf('(i) Naive continuation value == exponential (beta0=%g), should be zero: %2.8f \n',Params.beta0,max(abs(VnAalt(:)-Vexp(:))))

% (ii) beta0=1, Naive V and Valt equal exponential; (iii) beta0=1, Sophisticated V and Valt equal exponential
beta0_store=Params.beta0; Params.beta0=1;
[Vexp1,~]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsE);
[VnA1,~,VnA1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsN);
vfoptionsS=vfoptions; vfoptionsS.quasi_hyperbolic='Sophisticated';
[VsS1,~,VsS1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsS);
fprintf('(ii) Naive V == exponential (beta0=1), should be zero: %2.8f \n',max(abs(VnA1(:)-Vexp1(:))))
fprintf('(ii) Naive Valt == exponential (beta0=1), should be zero: %2.8f \n',max(abs(VnA1alt(:)-Vexp1(:))))
fprintf('(iii) Sophisticated V == exponential (beta0=1), should be zero: %2.8f \n',max(abs(VsS1(:)-Vexp1(:))))
fprintf('(iii) Sophisticated Valt == exponential (beta0=1), should be zero: %2.8f \n',max(abs(VsS1alt(:)-Vexp1(:))))
Params.beta0=beta0_store;

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
