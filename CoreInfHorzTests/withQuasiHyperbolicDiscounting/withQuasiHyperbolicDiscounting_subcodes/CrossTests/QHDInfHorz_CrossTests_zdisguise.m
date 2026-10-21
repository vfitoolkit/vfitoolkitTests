function output=QHDInfHorz_CrossTests_zdisguise(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-tests: run the same model disguised as a different z process, and check the answer follows.
%
% A: SPLIT. Replace each z state by two identical copies, and split the transition probabilities
%    evenly between them. That is the very same stochastic process, so V and Policy must simply
%    be duplicated column by column. This exercises the z loops and the pi_z indexing at a
%    different N_z, and it is the check that would catch a raw that indexed pi_z by the wrong
%    dimension, or summed over z rather than z'.
%
% B: RELABEL. Permute the z states (grid and transition matrix together). V and Policy must
%    permute with them. A raw that mixed up z and zprime anywhere passes A but fails B, because
%    a permuted pi_z is no longer symmetric in the way the split one is.
%
% Both are run for Naive and Sophisticated, with and without d, on the pure-discretization path.
% Sophisticated is a roughly-equal check throughout: it is pure VFI capped at maxiter=1000, and
% the split model has twice as many states, so the two runs stop at slightly different residuals.

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn_nod=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
ReturnFn_d=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);

vfoptionsQH=struct();
vfoptionsQH.exoticpreferences='QuasiHyperbolic';
vfoptionsQH.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;

% A: each z state split into two identical copies, transitions split evenly between the copies
n_z_split=2*N_z;
z_grid_split=repelem(z_grid,2,1);
pi_z_split=repelem(pi_z,2,2)/2;
fprintf('\n===== z-disguise cross-tests ===== \n')
fprintf('split pi_z rows sum to one, this should be zero: %2.8e \n',max(abs(sum(pi_z_split,2)-1)))

% B: relabelling, interleave-split so the permuted grid is clearly not monotone
perm=[(2:2:N_z)';(1:2:N_z)'];
z_grid_perm=z_grid(perm);
pi_z_perm=pi_z(perm,perm);
fprintf('permuted pi_z rows sum to one, this should be zero: %2.8e \n',max(abs(sum(pi_z_perm,2)-1)))

%% ---- without d, Naive
fprintf('--- nod, Naive --- \n')
vfN=vfoptionsQH; vfN.quasi_hyperbolic='Naive';
[V0,P0,V0alt,P0alt]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,[],vfN);
[Vsp,Psp,Vspalt,Pspalt]=ValueFnIter_InfHorz(0,n_a,n_z_split,[],a_grid,z_grid_split,pi_z_split,ReturnFn_nod,Params,DiscountFactorParamNames,[],vfN);
fprintf('A nod Naive: split V vs duplicated base V, this should be zero: %2.8e \n',max(abs(Vsp(:)-reshape(repelem(V0,1,2),[N_a*n_z_split,1]))))
fprintf('A nod Naive: split Policy vs duplicated base Policy, this should be zero: %.3e \n',max(abs(Psp(:)-reshape(repelem(P0,1,1,2),[N_a*n_z_split,1]))))
fprintf('A nod Naive: split Valt vs duplicated base Valt, this should be zero: %2.8e \n',max(abs(Vspalt(:)-reshape(repelem(V0alt,1,2),[N_a*n_z_split,1]))))
fprintf('A nod Naive: split Policyalt vs duplicated base Policyalt, this should be zero: %.3e \n',max(abs(Pspalt(:)-reshape(repelem(P0alt,1,1,2),[N_a*n_z_split,1]))))

[Vpm,Ppm,Vpmalt,Ppmalt]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid_perm,pi_z_perm,ReturnFn_nod,Params,DiscountFactorParamNames,[],vfN);
fprintf('B nod Naive: permuted V vs permuted base V, this should be zero: %2.8e \n',max(abs(Vpm(:)-reshape(V0(:,perm),[N_a*N_z,1]))))
fprintf('B nod Naive: permuted Policy vs permuted base Policy, this should be zero: %.3e \n',max(abs(Ppm(:)-reshape(P0(:,:,perm),[N_a*N_z,1]))))
fprintf('B nod Naive: permuted Valt vs permuted base Valt, this should be zero: %2.8e \n',max(abs(Vpmalt(:)-reshape(V0alt(:,perm),[N_a*N_z,1]))))
fprintf('B nod Naive: permuted Policyalt vs permuted base Policyalt, this should be zero: %.3e \n',max(abs(Ppmalt(:)-reshape(P0alt(:,:,perm),[N_a*N_z,1]))))

clear V0 P0 V0alt P0alt Vsp Psp Vspalt Pspalt Vpm Ppm Vpmalt Ppmalt

%% ---- without d, Sophisticated
fprintf('--- nod, Sophisticated --- \n')
vfS=vfoptionsQH; vfS.quasi_hyperbolic='Sophisticated';
[V0,P0,V0alt]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,[],vfS);
[Vsp,Psp,Vspalt]=ValueFnIter_InfHorz(0,n_a,n_z_split,[],a_grid,z_grid_split,pi_z_split,ReturnFn_nod,Params,DiscountFactorParamNames,[],vfS);
fprintf('A nod Sophisticated: split Vhat vs duplicated base Vhat, this should be near zero: %2.8e \n',max(abs(Vsp(:)-reshape(repelem(V0,1,2),[N_a*n_z_split,1]))))
fprintf('A nod Sophisticated: split Vunderbar vs duplicated base Vunderbar, this should be near zero: %2.8e \n',max(abs(Vspalt(:)-reshape(repelem(V0alt,1,2),[N_a*n_z_split,1]))))
fprintf('A nod Sophisticated: split Policy vs duplicated base Policy, this should be zero: %.3e \n',max(abs(Psp(:)-reshape(repelem(P0,1,1,2),[N_a*n_z_split,1]))))

[Vpm,Ppm,Vpmalt]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid_perm,pi_z_perm,ReturnFn_nod,Params,DiscountFactorParamNames,[],vfS);
fprintf('B nod Sophisticated: permuted Vhat vs permuted base Vhat, this should be near zero: %2.8e \n',max(abs(Vpm(:)-reshape(V0(:,perm),[N_a*N_z,1]))))
fprintf('B nod Sophisticated: permuted Vunderbar vs permuted base Vunderbar, this should be near zero: %2.8e \n',max(abs(Vpmalt(:)-reshape(V0alt(:,perm),[N_a*N_z,1]))))
fprintf('B nod Sophisticated: permuted Policy vs permuted base Policy, this should be zero: %.3e \n',max(abs(Ppm(:)-reshape(P0(:,:,perm),[N_a*N_z,1]))))

clear V0 P0 V0alt Vsp Psp Vspalt Vpm Ppm Vpmalt

%% ---- with d, Naive (refinement, the default)
fprintf('--- d, Naive --- \n')
vfN=vfoptionsQH; vfN.quasi_hyperbolic='Naive';
[V0,P0,V0alt,P0alt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,[],vfN);
[Vsp,Psp,Vspalt,Pspalt]=ValueFnIter_InfHorz(n_d,n_a,n_z_split,d_grid,a_grid,z_grid_split,pi_z_split,ReturnFn_d,Params,DiscountFactorParamNames,[],vfN);
fprintf('A d Naive: split V vs duplicated base V, this should be zero: %2.8e \n',max(abs(Vsp(:)-reshape(repelem(V0,1,2),[N_a*n_z_split,1]))))
fprintf('A d Naive: split Policy vs duplicated base Policy, this should be zero: %.3e \n',max(abs(Psp(:)-reshape(repelem(P0,1,1,2),[2*N_a*n_z_split,1]))))
fprintf('A d Naive: split Valt vs duplicated base Valt, this should be zero: %2.8e \n',max(abs(Vspalt(:)-reshape(repelem(V0alt,1,2),[N_a*n_z_split,1]))))
fprintf('A d Naive: split Policyalt vs duplicated base Policyalt, this should be zero: %.3e \n',max(abs(Pspalt(:)-reshape(repelem(P0alt,1,1,2),[2*N_a*n_z_split,1]))))

[Vpm,Ppm,Vpmalt,Ppmalt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid_perm,pi_z_perm,ReturnFn_d,Params,DiscountFactorParamNames,[],vfN);
fprintf('B d Naive: permuted V vs permuted base V, this should be zero: %2.8e \n',max(abs(Vpm(:)-reshape(V0(:,perm),[N_a*N_z,1]))))
fprintf('B d Naive: permuted Policy vs permuted base Policy, this should be zero: %.3e \n',max(abs(Ppm(:)-reshape(P0(:,:,perm),[2*N_a*N_z,1]))))
fprintf('B d Naive: permuted Valt vs permuted base Valt, this should be zero: %2.8e \n',max(abs(Vpmalt(:)-reshape(V0alt(:,perm),[N_a*N_z,1]))))
fprintf('B d Naive: permuted Policyalt vs permuted base Policyalt, this should be zero: %.3e \n',max(abs(Ppmalt(:)-reshape(P0alt(:,:,perm),[2*N_a*N_z,1]))))

clear V0 P0 V0alt P0alt Vsp Psp Vspalt Pspalt Vpm Ppm Vpmalt Ppmalt

%% ---- with d, Sophisticated (refinement, the default)
fprintf('--- d, Sophisticated --- \n')
vfS=vfoptionsQH; vfS.quasi_hyperbolic='Sophisticated';
[V0,P0,V0alt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,[],vfS);
[Vsp,Psp,Vspalt]=ValueFnIter_InfHorz(n_d,n_a,n_z_split,d_grid,a_grid,z_grid_split,pi_z_split,ReturnFn_d,Params,DiscountFactorParamNames,[],vfS);
fprintf('A d Sophisticated: split Vhat vs duplicated base Vhat, this should be near zero: %2.8e \n',max(abs(Vsp(:)-reshape(repelem(V0,1,2),[N_a*n_z_split,1]))))
fprintf('A d Sophisticated: split Vunderbar vs duplicated base Vunderbar, this should be near zero: %2.8e \n',max(abs(Vspalt(:)-reshape(repelem(V0alt,1,2),[N_a*n_z_split,1]))))
fprintf('A d Sophisticated: split Policy vs duplicated base Policy, this should be zero: %.3e \n',max(abs(Psp(:)-reshape(repelem(P0,1,1,2),[2*N_a*n_z_split,1]))))

[Vpm,Ppm,Vpmalt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid_perm,pi_z_perm,ReturnFn_d,Params,DiscountFactorParamNames,[],vfS);
fprintf('B d Sophisticated: permuted Vhat vs permuted base Vhat, this should be near zero: %2.8e \n',max(abs(Vpm(:)-reshape(V0(:,perm),[N_a*N_z,1]))))
fprintf('B d Sophisticated: permuted Vunderbar vs permuted base Vunderbar, this should be near zero: %2.8e \n',max(abs(Vpmalt(:)-reshape(V0alt(:,perm),[N_a*N_z,1]))))
fprintf('B d Sophisticated: permuted Policy vs permuted base Policy, this should be zero: %.3e \n',max(abs(Ppm(:)-reshape(P0(:,:,perm),[2*N_a*N_z,1]))))

%%
output=struct();

end
