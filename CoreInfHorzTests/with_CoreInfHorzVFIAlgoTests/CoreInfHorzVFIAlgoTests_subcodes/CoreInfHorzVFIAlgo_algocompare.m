function output=CoreInfHorzVFIAlgo_algocompare(label,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions)
% Shared engine for the InfHorz VFI-algorithm tests.
%
% Compares the different ValueFnIter_InfHorz algorithm options for one model.
% The value function iteration converges to the same fixed point regardless of
% the algorithm used, so all the variants below should give the same V and
% Policy (up to the convergence tolerance for V; Policy indices should match
% exactly). With the grid interpolation layer (GI) the same is true, and preGI
% vs postGI should also agree.
%
% Inputs: n_d/n_z (and d_grid/z_grid/pi_z) already zeroed as needed for this
% model; vfoptions already contains n_e/e_grid/pi_e for the with-e models, and
% is otherwise empty (no GI / howards fields set).
%
% Algorithm options being swept:
%   vfoptions.howardsgreedy   = 0 (iterated/modified-policy-iteration), 1 (greedy/policy-iteration),
%                               and (GI only) 2/3 (HowardMix: greedy on one grid, iterated on the other)
%   vfoptions.howardssparse   = 0 (indexed Howards) or 1 (sparse-matrix Howards)
%   vfoptions.howards         = number of Howards improvement steps (0 => pure VFI, no acceleration)
%   vfoptions.lowmemory       = 0 (precompute return matrix) or 1 (recompute, saves memory)
%   vfoptions.gridinterplayer = 0 or 1, and when 1: vfoptions.preGI = 0 (postGI) or 1 (preGI)
%
% NOTE on commented-out lines: some (option x GI) combinations are not yet
% implemented in the toolkit and would error. Per request they are included but
% commented out, so the intended coverage is documented in one place.

fprintf('\n================ %s ================\n',label);
DF=DiscountFactorParamNames;

%% Part 1: WITHOUT grid interpolation, all algorithm variants should give the same V and Policy
% Reference: howardsgreedy=0, howardssparse=0, lowmemory=0 (default howards acceleration)
vfo=vfoptions; vfo.howardsgreedy=0; vfo.howardssparse=0;
[V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);

% howardsgreedy=0, howardssparse=1
vfo=vfoptions; vfo.howardsgreedy=0; vfo.howardssparse=1;
[V2,Policy2]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('noGI: greedy0/sparse1 vs ref,          V   should be ~0: %2.8f \n',max(abs(V(:)-V2(:))));
fprintf('noGI: greedy0/sparse1 vs ref,          Pol should be  0: %2.8f \n',max(abs(Policy(:)-Policy2(:))));

% howardsgreedy=1, howardssparse=0
vfo=vfoptions; vfo.howardsgreedy=1; vfo.howardssparse=0;
[V3,Policy3]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('noGI: greedy1/sparse0 vs ref,          V   should be ~0: %2.8f \n',max(abs(V(:)-V3(:))));
fprintf('noGI: greedy1/sparse0 vs ref,          Pol should be  0: %2.8f \n',max(abs(Policy(:)-Policy3(:))));

% howardsgreedy=1, howardssparse=1 (greedy ignores sparse; included for completeness of the 2x2)
vfo=vfoptions; vfo.howardsgreedy=1; vfo.howardssparse=1;
[V4,Policy4]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('noGI: greedy1/sparse1 vs ref,          V   should be ~0: %2.8f \n',max(abs(V(:)-V4(:))));
fprintf('noGI: greedy1/sparse1 vs ref,          Pol should be  0: %2.8f \n',max(abs(Policy(:)-Policy4(:))));

% Howards OFF (howards=0 => pure value function iteration, no Howards acceleration)
vfo=vfoptions; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.howards=0;
[V5,Policy5]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('noGI: howards=0 (pure VFI) vs ref,     V   should be ~0: %2.8f \n',max(abs(V(:)-V5(:))));
fprintf('noGI: howards=0 (pure VFI) vs ref,     Pol should be  0: %2.8f \n',max(abs(Policy(:)-Policy5(:))));

% lowmemory=1 (recompute return matrix each iteration; loops over greedy=0 iterated Howards)
vfo=vfoptions; vfo.lowmemory=1; vfo.howardsgreedy=0; vfo.howardssparse=0;
[V6,Policy6]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('noGI: lowmemory1/sparse0 vs ref,       V   should be ~0: %2.8f \n',max(abs(V(:)-V6(:))));
fprintf('noGI: lowmemory1/sparse0 vs ref,       Pol should be  0: %2.8f \n',max(abs(Policy(:)-Policy6(:))));

% lowmemory=1, howardssparse=1 (sparse Howards is applied for the nod models; for the d models the
% lowmemory path ignores sparse, so this just repeats lowmemory1/sparse0 there -- harmless)
vfo=vfoptions; vfo.lowmemory=1; vfo.howardsgreedy=0; vfo.howardssparse=1;
[V7,Policy7]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('noGI: lowmemory1/sparse1 vs ref,       V   should be ~0: %2.8f \n',max(abs(V(:)-V7(:))));
fprintf('noGI: lowmemory1/sparse1 vs ref,       Pol should be  0: %2.8f \n',max(abs(Policy(:)-Policy7(:))));

clear V2 V3 V4 V5 V6 V7 Policy2 Policy3 Policy4 Policy5 Policy6 Policy7

%% Part 2: WITH grid interpolation, all (implemented) algorithm variants should give the same V and Policy
% Reference GI: howardsgreedy=0, howardssparse=0, preGI=0 (postGI, the default)
vfoGI=vfoptions; vfoGI.gridinterplayer=1; vfoGI.ngridinterp=5;
vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0;
[VG,PolicyG]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);

% howardsgreedy=1, howardssparse=0
vfo=vfoGI; vfo.howardsgreedy=1; vfo.howardssparse=0;
[VG3,PolicyG3]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   greedy1/sparse0 vs ref,          V   should be ~0: %2.8f \n',max(abs(VG(:)-VG3(:))));
fprintf('GI:   greedy1/sparse0 vs ref,          Pol should be  0: %2.8f \n',max(abs(PolicyG(:)-PolicyG3(:))));

% howardsgreedy=2 (HowardMix: greedy on a_grid, iterated on aprime_grid)
vfo=vfoGI; vfo.howardsgreedy=2; vfo.howardssparse=0;
[VG4,PolicyG4]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   greedy2 (HowardMix) vs ref,      V   should be ~0: %2.8f \n',max(abs(VG(:)-VG4(:))));
fprintf('GI:   greedy2 (HowardMix) vs ref,      Pol should be  0: %2.8f \n',max(abs(PolicyG(:)-PolicyG4(:))));

% howardsgreedy=3 (HowardMix2: iterated on a_grid, greedy on aprime_grid)
vfo=vfoGI; vfo.howardsgreedy=3; vfo.howardssparse=0;
[VG5,PolicyG5]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   greedy3 (HowardMix2) vs ref,     V   should be ~0: %2.8f \n',max(abs(VG(:)-VG5(:))));
fprintf('GI:   greedy3 (HowardMix2) vs ref,     Pol should be  0: %2.8f \n',max(abs(PolicyG(:)-PolicyG5(:))));

% howardssparse=1 with GI: implemented only for lowmemory=1 (and only postGI)
vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=1; vfo.lowmemory=1;
[VG6,PolicyG6]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   sparse1/lowmemory1 vs ref,       V   should be ~0: %2.8f \n',max(abs(VG(:)-VG6(:))));
fprintf('GI:   sparse1/lowmemory1 vs ref,       Pol should be  0: %2.8f \n',max(abs(PolicyG(:)-PolicyG6(:))));

% howardsgreedy=0, howardssparse=1, lowmemory=0 with GI -- NOT YET IMPLEMENTED (only lowmemory=1).
% (ValueFnIter_InfHorz_GridInterpLayer.m: 'howardssparse=1 only implemented for lowmemory=1'). Included, commented out:
% vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=1; vfo.lowmemory=0;
% [VG2,PolicyG2]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
% fprintf('GI:   greedy0/sparse1/lowmemory0 vs ref, V should be ~0: %2.8f \n',max(abs(VG(:)-VG2(:))));

clear VG3 VG4 VG5 VG6 PolicyG3 PolicyG4 PolicyG5 PolicyG6

%% Part 3: WITH grid interpolation, preGI and postGI should give the same V and Policy
% postGI (preGI=0) reference is VG (greedy0) from Part 2.
% preGI, greedy=0, sparse=0
vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.preGI=1;
[VpreA,PolicypreA]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   preGI vs postGI (greedy0),       V   should be ~0: %2.8f \n',max(abs(VG(:)-VpreA(:))));
fprintf('GI:   preGI vs postGI (greedy0),       Pol should be  0: %2.8f \n',max(abs(PolicyG(:)-PolicypreA(:))));

% preGI, greedy=1, sparse=0
vfo=vfoGI; vfo.howardsgreedy=1; vfo.howardssparse=0; vfo.preGI=1;
[VpreB,PolicypreB]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   preGI vs postGI (greedy1),       V   should be ~0: %2.8f \n',max(abs(VG(:)-VpreB(:))));
fprintf('GI:   preGI vs postGI (greedy1),       Pol should be  0: %2.8f \n',max(abs(PolicyG(:)-PolicypreB(:))));

% preGI, greedy=2 (HowardMix) -- preGI implements greedy=2
vfo=vfoGI; vfo.howardsgreedy=2; vfo.howardssparse=0; vfo.preGI=1;
[VpreC,PolicypreC]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   preGI vs postGI (greedy2),       V   should be ~0: %2.8f \n',max(abs(VG(:)-VpreC(:))));
fprintf('GI:   preGI vs postGI (greedy2),       Pol should be  0: %2.8f \n',max(abs(PolicyG(:)-PolicypreC(:))));

% preGI, greedy=3 -- NOT YET IMPLEMENTED for preGI (postGI only). Included, commented out:
% vfo=vfoGI; vfo.howardsgreedy=3; vfo.howardssparse=0; vfo.preGI=1;
% [VpreD,PolicypreD]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
% fprintf('GI:   preGI vs postGI (greedy3),       V   should be ~0: %2.8f \n',max(abs(VG(:)-VpreD(:))));

% preGI, sparse=1 -- NOT YET IMPLEMENTED for preGI (postGI+lowmemory=1 only). Included, commented out:
% vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=1; vfo.lowmemory=1; vfo.preGI=1;
% [VpreE,PolicypreE]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
% fprintf('GI:   preGI vs postGI (sparse1),       V   should be ~0: %2.8f \n',max(abs(VG(:)-VpreE(:))));

clear VG PolicyG VpreA VpreB VpreC PolicypreA PolicypreB PolicypreC

%% Part 4: big n_a (=1500), with and without grid interpolation should give very similar V
% (Not exactly equal: GI refines the choice between grid points, so V_GI is a
%  touch higher; with 1500 points the difference is tiny.)
vfo=vfoptions; % without GI
[Vbig,~]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
vfo=vfoptions; vfo.gridinterplayer=1; vfo.ngridinterp=5; % with GI
[VbigGI,~]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('big n_a=%d: with vs without GI, V should be very similar (small): %2.8f \n',n_a_big,max(abs(Vbig(:)-VbigGI(:))));

%% Part 5 (DIAGNOSTIC): is the postGI-vs-preGI gap caused by postGI's search window?
% preGI searches the ENTIRE fine aprime grid (the exact optimizer of the GI objective).
% postGI only searches +-maxaprimediff COARSE cells around the coarse-grid optimum. If that
% window is the cause of the gap, widening maxaprimediff should make postGI converge to preGI.
% If it does NOT converge, the gap is a genuine bug rather than the search window.
% (Note: a fine-grid no-GI reference does NOT work here -- GI keeps V on the coarse grid and
%  linearly interpolates the continuation, so BOTH preGI and postGI share an irreducible
%  O(coarse-grid) interpolation error and sit ~equally far from any fine no-GI solution.)
vfoGI=vfoptions; vfoGI.gridinterplayer=1; vfoGI.ngridinterp=5;
vfo=vfoGI; vfo.preGI=1;                       [Vpre,~]  =ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
vfo=vfoGI; vfo.preGI=0;                        [Vpost,~] =ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); % default window
fprintf('postGI default window: |preGI-postGI| = %2.8f,  min(Vpre-Vpost) = %2.8f  (>=0 => preGI weakly higher, i.e. postGI under-searches) \n',max(abs(Vpre(:)-Vpost(:))),min(Vpre(:)-Vpost(:)));
mad=floor((n_a-1)/2); % widest VALID window: postGI clamps aprimeshifter to [1+mad, N_a-mad], so this covers the whole grid = full search
vfo=vfoGI; vfo.preGI=0; vfo.maxaprimediff=mad; [VpostW,~] =ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('postGI WIDE window (maxaprimediff=%d): |preGI-postGI_wide| = %2.8f  (if ~0 => the window was the cause; preGI is the exact reference) \n',mad,max(abs(Vpre(:)-VpostW(:))));

output=struct();

end
