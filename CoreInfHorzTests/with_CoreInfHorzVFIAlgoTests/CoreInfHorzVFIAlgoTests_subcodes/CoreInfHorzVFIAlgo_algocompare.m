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
% commented out, so the intended coverage is documented in one place. Two such combos remain:
% preGI with howardsgreedy=3, and preGI with howardssparse=1 (both are postGI-only).

fprintf('\n================ %s ================\n',label);
DF=DiscountFactorParamNames;

%% Part 1: WITHOUT grid interpolation, all algorithm variants should give the same V and Policy
% Reference: howardsgreedy=0, howardssparse=0, lowmemory=0 (default howards acceleration)
vfo=vfoptions; vfo.howardsgreedy=0; vfo.howardssparse=0;
[V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);

% howardsgreedy=0, howardssparse=1
vfo=vfoptions; vfo.howardsgreedy=0; vfo.howardssparse=1;
[V2,Policy2]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('noGI: greedy0/sparse1 vs ref,          V   should be ~0: %.3e \n',max(abs(V(:)-V2(:))));
fprintf('noGI: greedy0/sparse1 vs ref,          Pol should be  0: %.3e \n',max(abs(Policy(:)-Policy2(:))));

% howardsgreedy=1, howardssparse=0
vfo=vfoptions; vfo.howardsgreedy=1; vfo.howardssparse=0;
[V3,Policy3]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('noGI: greedy1/sparse0 vs ref,          V   should be ~0: %.3e \n',max(abs(V(:)-V3(:))));
fprintf('noGI: greedy1/sparse0 vs ref,          Pol should be  0: %.3e \n',max(abs(Policy(:)-Policy3(:))));

% howardsgreedy=1, howardssparse=1 (greedy ignores sparse; included for completeness of the 2x2)
vfo=vfoptions; vfo.howardsgreedy=1; vfo.howardssparse=1;
[V4,Policy4]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('noGI: greedy1/sparse1 vs ref,          V   should be ~0: %.3e \n',max(abs(V(:)-V4(:))));
fprintf('noGI: greedy1/sparse1 vs ref,          Pol should be  0: %.3e \n',max(abs(Policy(:)-Policy4(:))));

% Howards OFF (howards=0 => pure value function iteration, no Howards acceleration)
vfo=vfoptions; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.howards=0;
[V5,Policy5]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('noGI: howards=0 (pure VFI) vs ref,     V   should be ~0: %.3e \n',max(abs(V(:)-V5(:))));
fprintf('noGI: howards=0 (pure VFI) vs ref,     Pol should be  0: %.3e \n',max(abs(Policy(:)-Policy5(:))));

% lowmemory=1. In infinite horizon this exists only for refinement, which needs a decision variable
% d: refinement builds the refined return matrix one z at a time. A model with no d has no lowmemory
% option at all (ValueFnIter_InfHorz_PureDiscretization errors if you ask for one), so skip these two.
if prod(n_d)>0
    vfo=vfoptions; vfo.lowmemory=1; vfo.howardsgreedy=0; vfo.howardssparse=0;
    [V6,Policy6]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    fprintf('noGI: lowmemory1/sparse0 vs ref,       V   should be ~0: %.3e \n',max(abs(V(:)-V6(:))));
    fprintf('noGI: lowmemory1/sparse0 vs ref,       Pol should be  0: %.3e \n',max(abs(Policy(:)-Policy6(:))));

    % lowmemory=1, howardssparse=1. Refinement hands the refined (nod-shaped) matrix to the same raws
    % as lowmemory=0, so howardssparse still applies here and this is a genuinely different solve.
    vfo=vfoptions; vfo.lowmemory=1; vfo.howardsgreedy=0; vfo.howardssparse=1;
    [V7,Policy7]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    fprintf('noGI: lowmemory1/sparse1 vs ref,       V   should be ~0: %.3e \n',max(abs(V(:)-V7(:))));
    fprintf('noGI: lowmemory1/sparse1 vs ref,       Pol should be  0: %.3e \n',max(abs(Policy(:)-Policy7(:))));
else
    fprintf('noGI: lowmemory1 not run, this model has no d so infinite horizon has no lowmemory option \n');
end

% solnmethod='purediscretization'. With a d the default is refinement, which takes the max over d first and
% then hands a nod-shaped return matrix to the nod raws, so the with-d raws are only reached this way.
if prod(n_d)>0
    vfo=vfoptions; vfo.solnmethod='purediscretization';
    [V8,Policy8]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    fprintf('noGI: purediscretization vs ref,       V   should be ~0: %.3e \n',max(abs(V(:)-V8(:))));
    fprintf('noGI: purediscretization vs ref,       Pol should be  0: %.3e \n',max(abs(Policy(:)-Policy8(:))));

    % and with howards=0, so the with-d raws' own Howards loop is checked against pure VFI
    vfo=vfoptions; vfo.solnmethod='purediscretization'; vfo.howards=0;
    [V9,Policy9]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    fprintf('noGI: purediscretization howards=0,    V   should be ~0: %.3e \n',max(abs(V(:)-V9(:))));
    fprintf('noGI: purediscretization howards=0,    Pol should be  0: %.3e \n',max(abs(Policy(:)-Policy9(:))));
end

% outputkron=1. With one d and one a the un-Kron'd Policy is the Kron one with its dimensions spelt out, so the two
% must hold the same numbers. With a d this is the 2-row Policy, which the outputkron=1 reshape has to keep both rows of.
% (Only the models without e and without semiz, those are the outputkron branches checked here.)
if ~isfield(vfoptions,'n_e') && ~isfield(vfoptions,'n_semiz')
    vfo=vfoptions; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.outputkron=1;
    [V10,Policy10]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    fprintf('noGI: outputkron=1 vs ref,             V   should be  0: %.3e \n',max(abs(V(:)-V10(:))));
    fprintf('noGI: outputkron=1 vs ref,             Pol should be  0: %.3e \n',max(abs(Policy(:)-Policy10(:))));
    if prod(n_d)>0
        vfo=vfoptions; vfo.solnmethod='purediscretization'; vfo.outputkron=1;
        [V11,Policy11]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
        fprintf('noGI: purediscretization outputkron=1, V   should be ~0: %.3e \n',max(abs(V(:)-V11(:))));
        fprintf('noGI: purediscretization outputkron=1, Pol should be  0: %.3e \n',max(abs(Policy(:)-Policy11(:))));
    end
end

% endogenousvfi=1 (endogenous VFI; relative VFI when there is no z). Iterates exactly as VFI does but stops on the span over a
% of V-Vold within each exogenous state, then recovers V in levels, so V agrees with the reference to the same order as the
% reference is from the true V (beta/(1-beta)*tolerance each), and Policy exactly. Each solver gets its own endogenous twin.
% With semiz the span is over (a,semiz) within each (z,e), as the semiz transitions depend on d.
vfo=vfoptions; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.endogenousvfi=1;
[VE1,PolicyE1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('noGI: endoVFI greedy0/sparse0 vs ref,  V   should be ~0: %.3e \n',max(abs(V(:)-VE1(:))));
fprintf('noGI: endoVFI greedy0/sparse0 vs ref,  Pol should be  0: %.3e \n',max(abs(Policy(:)-PolicyE1(:))));

vfo=vfoptions; vfo.howardsgreedy=0; vfo.howardssparse=1; vfo.endogenousvfi=1;
[VE2,PolicyE2]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('noGI: endoVFI greedy0/sparse1 vs ref,  V   should be ~0: %.3e \n',max(abs(V(:)-VE2(:))));
fprintf('noGI: endoVFI greedy0/sparse1 vs ref,  Pol should be  0: %.3e \n',max(abs(Policy(:)-PolicyE2(:))));

vfo=vfoptions; vfo.howardsgreedy=1; vfo.howardssparse=0; vfo.endogenousvfi=1;
[VE3,PolicyE3]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('noGI: endoVFI greedy1 vs ref,          V   should be ~0: %.3e \n',max(abs(V(:)-VE3(:))));
fprintf('noGI: endoVFI greedy1 vs ref,          Pol should be  0: %.3e \n',max(abs(Policy(:)-PolicyE3(:))));

vfo=vfoptions; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.howards=0; vfo.endogenousvfi=1;
[VE4,PolicyE4]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('noGI: endoVFI howards=0 vs ref,        V   should be ~0: %.3e \n',max(abs(V(:)-VE4(:))));
fprintf('noGI: endoVFI howards=0 vs ref,        Pol should be  0: %.3e \n',max(abs(Policy(:)-PolicyE4(:))));

if prod(n_d)>0
    vfo=vfoptions; vfo.lowmemory=1; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.endogenousvfi=1;
    [VE5,PolicyE5]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    fprintf('noGI: endoVFI lowmemory1 vs ref,       V   should be ~0: %.3e \n',max(abs(V(:)-VE5(:))));
    fprintf('noGI: endoVFI lowmemory1 vs ref,       Pol should be  0: %.3e \n',max(abs(Policy(:)-PolicyE5(:))));

    % the with-d raws are only reached with solnmethod='purediscretization'
    vfo=vfoptions; vfo.solnmethod='purediscretization'; vfo.endogenousvfi=1;
    [VE6,PolicyE6]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    fprintf('noGI: endoVFI purediscretization,      V   should be ~0: %.3e \n',max(abs(V(:)-VE6(:))));
    fprintf('noGI: endoVFI purediscretization,      Pol should be  0: %.3e \n',max(abs(Policy(:)-PolicyE6(:))));

    vfo=vfoptions; vfo.solnmethod='purediscretization'; vfo.howards=0; vfo.endogenousvfi=1;
    [VE7,PolicyE7]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    fprintf('noGI: endoVFI purediscretization h=0,  V   should be ~0: %.3e \n',max(abs(V(:)-VE7(:))));
    fprintf('noGI: endoVFI purediscretization h=0,  Pol should be  0: %.3e \n',max(abs(Policy(:)-PolicyE7(:))));
end

% The discriminator: without it, a twin that was just plain VFI would pass every check above. Start from the endogenous
% solution plus an arbitrary function of the exogenous states, f, which changes neither the policy nor the change in V
% net of a function of the exogenous states. The residual of VE1 has a span within each exogenous state of at most
% beta*tolerance, so endogenous VFI must stop after the first iteration, and recover the same V in levels. Plain VFI from
% the same start is out by about (1-beta)*f in levels, so must run into maxiter=3.
% (V is [n_a,(n_semiz),(n_z),(n_e)]; f varies over the exogenous dimensions only, so not over semiz)
szV=size(V); if isfield(vfoptions,'n_semiz'), kexo=3; else, kexo=2; end
fexo=reshape(10+3*cos(1:prod(szV(kexo:end))),[ones(1,kexo-1),szV(kexo:end)]);
vfo=vfoptions; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.endogenousvfi=1; vfo.maxiter=3; vfo.V0=VE1+fexo;
lastwarn('');
[VE8,PolicyE8]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
msgE8=lastwarn;
fprintf('noGI: endoVFI from VE+f(exog), maxiter=3, converged,     should be one: %i \n',~contains(msgE8,'maximum number of iterations'));
fprintf('noGI: endoVFI from VE+f(exog), maxiter=3 vs ref,  V   should be ~0: %.3e \n',max(abs(V(:)-VE8(:))));
fprintf('noGI: endoVFI from VE+f(exog), maxiter=3 vs ref,  Pol should be  0: %.3e \n',max(abs(Policy(:)-PolicyE8(:))));
vfo=vfoptions; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.maxiter=3; vfo.V0=VE1+fexo; vfo.endogenousvfi=0;
lastwarn('');
[~,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
msgE9=lastwarn;
fprintf('noGI: plain VFI from VE+f(exog), maxiter=3, hit maxiter, should be one: %i \n',contains(msgE9,'maximum number of iterations'));

% The defaults that depend on the discount factor (see ValueFnIter_InfHorz): from beta=0.97 endogenous VFI with 40 Howards
% steps (indexed Howards, as at these grid sizes without GI), below it plain VFI with round(2/(1-beta)) Howards steps. So the
% default solve must be identical to the explicit one: the same code runs.
ParamsDefault=Params;
for bdef=[0.95,0.99]
    ParamsDefault.(DF{1})=bdef;
    vfo=vfoptions;
    [VD,PolicyD]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsDefault,DF,[],vfo);
    vfo=vfoptions; vfo.howardssparse=0;
    if bdef>=0.97
        vfo.endogenousvfi=1; vfo.howards=40;
    else
        vfo.endogenousvfi=0; vfo.howards=round(2/(1-bdef));
    end
    [VX,PolicyX]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsDefault,DF,[],vfo);
    fprintf('noGI: defaults at beta=%1.2f vs explicit (endogenousvfi=%i, howards=%i), V   should be  0: %.3e \n',bdef,vfo.endogenousvfi,vfo.howards,max(abs(VD(:)-VX(:))));
    fprintf('noGI: defaults at beta=%1.2f vs explicit (endogenousvfi=%i, howards=%i), Pol should be  0: %.3e \n',bdef,vfo.endogenousvfi,vfo.howards,max(abs(PolicyD(:)-PolicyX(:))));
end
clear ParamsDefault VD VX PolicyD PolicyX

clear V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 Policy2 Policy3 Policy4 Policy5 Policy6 Policy7 Policy8 Policy9 Policy10 Policy11
clear VE1 VE2 VE3 VE4 VE5 VE6 VE7 VE8 PolicyE1 PolicyE2 PolicyE3 PolicyE4 PolicyE5 PolicyE6 PolicyE7 PolicyE8

%% Without z and without e, the grid interpolation layer is a single postGI solve that reads neither howardsgreedy nor
% howardssparse, so the Part 2 comparisons would compare a solve with itself, and preGI (Part 3) is not implemented.
% So these models stop after Part 1.
if prod(n_z)==0 && ~isfield(vfoptions,'n_e')
    fprintf('GI parts (2-5) not run, without z and e the grid interpolation layer is a single postGI solve (no Howards options, no preGI) \n');
    output=struct();
    return
end

%% With semiz, the grid interpolation layer is so far only postGI, so the semiz models run their own short GI part here
% and stop: HowardMix (2/3), postGIrepeat and preGI are refused. Sparse and greedy Howards are only implemented without
% a d1 (with a d1 howardssparse is ignored, it would be the same indexed solve so comparing it against the reference
% would test nothing, and howardsgreedy=1 is refused), so the models with a d1 skip those two.
if isfield(vfoptions,'n_semiz')
    vfoGI=vfoptions; vfoGI.gridinterplayer=1; vfoGI.ngridinterp=5;
    vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0;
    [VG,PolicyG]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);

    % Howards OFF (howards=0 => pure value function iteration on both the coarse and the fine stage)
    vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.howards=0;
    [VG8,PolicyG8]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    fprintf('GI:   howards=0 (pure VFI) vs ref,     V   should be ~0: %.3e \n',max(abs(VG(:)-VG8(:))));
    fprintf('GI:   howards=0 (pure VFI) vs ref,     Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyG8(:))));

    if length(n_d)>1
        % lowmemory=1 (with a d1 the fine return matrix is built one exogenous state, or with e one e, at a time)
        vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.lowmemory=1;
        [VG9,PolicyG9]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
        fprintf('GI:   lowmemory1 vs ref,               V   should be ~0: %.3e \n',max(abs(VG(:)-VG9(:))));
        fprintf('GI:   lowmemory1 vs ref,               Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyG9(:))));
    else
        % howardssparse=1 (without a d1; also the default with the grid interpolation layer)
        vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=1;
        [VG6,PolicyG6]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
        fprintf('GI:   sparse1 vs ref,                  V   should be ~0: %.3e \n',max(abs(VG(:)-VG6(:))));
        fprintf('GI:   sparse1 vs ref,                  Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyG6(:))));

        % howardsgreedy=1 (without a d1)
        vfo=vfoGI; vfo.howardsgreedy=1; vfo.howardssparse=0;
        [VG3,PolicyG3]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
        fprintf('GI:   greedy1 vs ref,                  V   should be ~0: %.3e \n',max(abs(VG(:)-VG3(:))));
        fprintf('GI:   greedy1 vs ref,                  Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyG3(:))));
    end
    % endogenousvfi=1 with semiz and GI: the semiz GI endogenous twins, for each config above
    vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.endogenousvfi=1;
    [VGE,PolicyGE]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    fprintf('GI:   endoVFI vs ref,                  V   should be ~0: %.3e \n',max(abs(VG(:)-VGE(:))));
    fprintf('GI:   endoVFI vs ref,                  Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyGE(:))));
    vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.howards=0; vfo.endogenousvfi=1;
    [VGE,PolicyGE]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    fprintf('GI:   endoVFI howards=0 vs ref,        V   should be ~0: %.3e \n',max(abs(VG(:)-VGE(:))));
    fprintf('GI:   endoVFI howards=0 vs ref,        Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyGE(:))));
    if length(n_d)>1
        vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.lowmemory=1; vfo.endogenousvfi=1;
        [VGE,PolicyGE]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
        fprintf('GI:   endoVFI lowmemory1 vs ref,       V   should be ~0: %.3e \n',max(abs(VG(:)-VGE(:))));
        fprintf('GI:   endoVFI lowmemory1 vs ref,       Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyGE(:))));
    else
        vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=1; vfo.endogenousvfi=1;
        [VGE,PolicyGE]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
        fprintf('GI:   endoVFI sparse1 vs ref,          V   should be ~0: %.3e \n',max(abs(VG(:)-VGE(:))));
        fprintf('GI:   endoVFI sparse1 vs ref,          Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyGE(:))));
        vfo=vfoGI; vfo.howardsgreedy=1; vfo.howardssparse=0; vfo.endogenousvfi=1;
        [VGE,PolicyGE]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
        fprintf('GI:   endoVFI greedy1 vs ref,          V   should be ~0: %.3e \n',max(abs(VG(:)-VGE(:))));
        fprintf('GI:   endoVFI greedy1 vs ref,          Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyGE(:))));
    end
    clear VG VG3 VG6 VG8 VG9 VGE PolicyG PolicyG3 PolicyG6 PolicyG8 PolicyG9 PolicyGE
    fprintf('GI parts (2-5) otherwise not run, with semiz the grid interpolation layer is so far only postGI \n');
    output=struct();
    return
end

%% Part 2: WITH grid interpolation, all (implemented) algorithm variants should give the same V and Policy
% Reference GI: howardsgreedy=0, howardssparse=0, preGI=0 (postGI, the default)
vfoGI=vfoptions; vfoGI.gridinterplayer=1; vfoGI.ngridinterp=5;
vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0;
[VG,PolicyG]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);

% howardsgreedy=1, howardssparse=0
vfo=vfoGI; vfo.howardsgreedy=1; vfo.howardssparse=0;
[VG3,PolicyG3]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   greedy1/sparse0 vs ref,          V   should be ~0: %.3e \n',max(abs(VG(:)-VG3(:))));
fprintf('GI:   greedy1/sparse0 vs ref,          Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyG3(:))));

% With e, howardsgreedy=2 and 3 (HowardMix) are not yet implemented for the grid interpolation layer, so the e models skip these two
if ~isfield(vfoptions,'n_e')
% howardsgreedy=2 (HowardMix: greedy on a_grid, iterated on aprime_grid)
vfo=vfoGI; vfo.howardsgreedy=2; vfo.howardssparse=0;
[VG4,PolicyG4]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   greedy2 (HowardMix) vs ref,      V   should be ~0: %.3e \n',max(abs(VG(:)-VG4(:))));
fprintf('GI:   greedy2 (HowardMix) vs ref,      Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyG4(:))));

% howardsgreedy=3 (HowardMix2: iterated on a_grid, greedy on aprime_grid)
vfo=vfoGI; vfo.howardsgreedy=3; vfo.howardssparse=0;
[VG5,PolicyG5]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   greedy3 (HowardMix2) vs ref,     V   should be ~0: %.3e \n',max(abs(VG(:)-VG5(:))));
fprintf('GI:   greedy3 (HowardMix2) vs ref,     Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyG5(:))));
end % with e, howardsgreedy=2 and 3 are not yet implemented for the grid interpolation layer

% howardssparse=1 with GI (postGI only; preGI sparse is still not implemented)
vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=1; vfo.lowmemory=1;
[VG6,PolicyG6]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   sparse1/lowmemory1 vs ref,       V   should be ~0: %.3e \n',max(abs(VG(:)-VG6(:))));
fprintf('GI:   sparse1/lowmemory1 vs ref,       Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyG6(:))));

% Same again with lowmemory=0. For the d models lowmemory only changes how the two (refined) return
% matrices get built -- one z at a time, instead of whole -- so it must not change the answer. For
% the nod models lowmemory has no effect on this code path at all.
vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=1; vfo.lowmemory=0;
[VG2,PolicyG2]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   sparse1/lowmemory0 vs ref,       V   should be ~0: %.3e \n',max(abs(VG(:)-VG2(:))));
fprintf('GI:   sparse1/lowmemory0 vs ref,       Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyG2(:))));

% And the two directly against each other. Unlike the comparisons against the reference above, this
% one is not merely '~0': it is the same arithmetic in a different construction order, and from the
% first while-loop onwards it is literally the same code, so it should be exactly zero.
fprintf('GI:   sparse1 lowmem0 vs lowmem1,      V   should be   0: %.3e \n',max(abs(VG6(:)-VG2(:))));
fprintf('GI:   sparse1 lowmem0 vs lowmem1,      Pol should be   0: %.3e \n',max(abs(PolicyG6(:)-PolicyG2(:))));

% postGIrepeat=1: a second post-GI layer, with the +-maxaprimediff window re-centred on the first layer's solution. At
% the default maxaprimediff the first layer already reaches the fine-grid optimum (see the postGIrepeat scans in
% doPart(4)), so the extra layer must not change the answer. This is what runs the repeat-layer code at all.
vfo=vfoGI; vfo.postGIrepeat=1;
[VG7,PolicyG7]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   postGIrepeat=1 vs ref,           V   should be ~0: %.3e \n',max(abs(VG(:)-VG7(:))));
fprintf('GI:   postGIrepeat=1 vs ref,           Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyG7(:))));

% endogenousvfi=1 with GI: each GI solver has its own endogenous twin, which switches from the coarse to the fine grid on
% the span too, and recovers V in levels after each fine-grid loop (the postGIrepeat layers included).
vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.endogenousvfi=1;
[VGE1,PolicyGE1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   endoVFI greedy0/sparse0 vs ref,  V   should be ~0: %.3e \n',max(abs(VG(:)-VGE1(:))));
fprintf('GI:   endoVFI greedy0/sparse0 vs ref,  Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyGE1(:))));

vfo=vfoGI; vfo.howardsgreedy=1; vfo.howardssparse=0; vfo.endogenousvfi=1;
[VGE2,PolicyGE2]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   endoVFI greedy1 vs ref,          V   should be ~0: %.3e \n',max(abs(VG(:)-VGE2(:))));
fprintf('GI:   endoVFI greedy1 vs ref,          Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyGE2(:))));

if ~isfield(vfoptions,'n_e') % HowardMix is not implemented with e
    vfo=vfoGI; vfo.howardsgreedy=2; vfo.howardssparse=0; vfo.endogenousvfi=1;
    [VGE3,PolicyGE3]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    fprintf('GI:   endoVFI greedy2 (HowardMix),     V   should be ~0: %.3e \n',max(abs(VG(:)-VGE3(:))));
    fprintf('GI:   endoVFI greedy2 (HowardMix),     Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyGE3(:))));

    vfo=vfoGI; vfo.howardsgreedy=3; vfo.howardssparse=0; vfo.endogenousvfi=1;
    [VGE4,PolicyGE4]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    fprintf('GI:   endoVFI greedy3 (HowardMix2),    V   should be ~0: %.3e \n',max(abs(VG(:)-VGE4(:))));
    fprintf('GI:   endoVFI greedy3 (HowardMix2),    Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyGE4(:))));
end

vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=1; vfo.lowmemory=1; vfo.endogenousvfi=1;
[VGE5,PolicyGE5]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   endoVFI sparse1/lowmemory1,      V   should be ~0: %.3e \n',max(abs(VG(:)-VGE5(:))));
fprintf('GI:   endoVFI sparse1/lowmemory1,      Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyGE5(:))));

vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=1; vfo.lowmemory=0; vfo.endogenousvfi=1;
[VGE6,PolicyGE6]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   endoVFI sparse1/lowmemory0,      V   should be ~0: %.3e \n',max(abs(VG(:)-VGE6(:))));
fprintf('GI:   endoVFI sparse1/lowmemory0,      Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyGE6(:))));

vfo=vfoGI; vfo.postGIrepeat=1; vfo.endogenousvfi=1;
[VGE7,PolicyGE7]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   endoVFI postGIrepeat=1,          V   should be ~0: %.3e \n',max(abs(VG(:)-VGE7(:))));
fprintf('GI:   endoVFI postGIrepeat=1,          Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicyGE7(:))));

clear VG2 VG3 VG4 VG5 VG6 VG7 PolicyG2 PolicyG3 PolicyG4 PolicyG5 PolicyG6 PolicyG7
clear VGE1 VGE2 VGE3 VGE4 VGE5 VGE6 VGE7 PolicyGE1 PolicyGE2 PolicyGE3 PolicyGE4 PolicyGE5 PolicyGE6 PolicyGE7

%% With e, the grid interpolation layer is so far only implemented for postGI with one endogenous state, so the e
% models stop here: Parts 3-5 are preGI, two endogenous states and the big grid (remove this block once they are)
if isfield(vfoptions,'n_e')
    fprintf('GI parts 3-5 not run, with e the grid interpolation layer is so far only postGI with one endogenous state \n');
    output=struct();
    return
end

%% Part 3: WITH grid interpolation, preGI and postGI should give the same V and Policy
% postGI (preGI=0) reference is VG (greedy0) from Part 2.
% preGI, greedy=0, sparse=0
vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.preGI=1;
[VpreA,PolicypreA]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   preGI vs postGI (greedy0),       V   should be ~0: %.3e \n',max(abs(VG(:)-VpreA(:))));
fprintf('GI:   preGI vs postGI (greedy0),       Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicypreA(:))));

% preGI, greedy=1, sparse=0
vfo=vfoGI; vfo.howardsgreedy=1; vfo.howardssparse=0; vfo.preGI=1;
[VpreB,PolicypreB]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   preGI vs postGI (greedy1),       V   should be ~0: %.3e \n',max(abs(VG(:)-VpreB(:))));
fprintf('GI:   preGI vs postGI (greedy1),       Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicypreB(:))));

% preGI, greedy=2 (HowardMix) -- preGI implements greedy=2
vfo=vfoGI; vfo.howardsgreedy=2; vfo.howardssparse=0; vfo.preGI=1;
[VpreC,PolicypreC]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('GI:   preGI vs postGI (greedy2),       V   should be ~0: %.3e \n',max(abs(VG(:)-VpreC(:))));
fprintf('GI:   preGI vs postGI (greedy2),       Pol should be  0: %.3e \n',max(abs(PolicyG(:)-PolicypreC(:))));

% endogenousvfi=1 with preGI (the preGI endogenous twins), against the postGI reference
for gg=0:2
    vfo=vfoGI; vfo.howardsgreedy=gg; vfo.howardssparse=0; vfo.preGI=1; vfo.endogenousvfi=1;
    [VpreE,PolicypreE]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    fprintf('GI:   endoVFI preGI vs postGI (greedy%i), V   should be ~0: %.3e \n',gg,max(abs(VG(:)-VpreE(:))));
    fprintf('GI:   endoVFI preGI vs postGI (greedy%i), Pol should be  0: %.3e \n',gg,max(abs(PolicyG(:)-PolicypreE(:))));
end
clear VpreE PolicypreE

% preGI, greedy=3 -- NOT YET IMPLEMENTED for preGI (postGI only). Included, commented out:
% vfo=vfoGI; vfo.howardsgreedy=3; vfo.howardssparse=0; vfo.preGI=1;
% [VpreD,PolicypreD]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
% fprintf('GI:   preGI vs postGI (greedy3),       V   should be ~0: %.3e \n',max(abs(VG(:)-VpreD(:))));

% preGI, sparse=1 -- NOT YET IMPLEMENTED for preGI (postGI+lowmemory=1 only). Included, commented out:
% vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=1; vfo.lowmemory=1; vfo.preGI=1;
% [VpreE,PolicypreE]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
% fprintf('GI:   preGI vs postGI (sparse1),       V   should be ~0: %.3e \n',max(abs(VG(:)-VpreE(:))));

%% Part 3b: the same preGI vs postGI check, but with TWO endogenous states
% This builds its own small model, because the one the caller passed in has a single endogenous
% state. It mirrors the caller's d-ness, so the 'nod' call covers the 2A nod raws and the 'd' call
% covers the 2A Refine raws, with neither repeating the other's work.
% Two things are being checked, each against a reference that is independent of it:
%   - preGI2A vs postGI2A. The preGI2A raws are the only two endogenous state grid interpolation
%     raws whose Howards code has not been touched, so this is the outside check on the postGI2A
%     ones (whose fine-grid stage had no working Howards at all until recently).
%   - postGI2A howardssparse=1 vs howardssparse=0. This is the check on the sparse raws, and in
%     particular on the joint index a1+N_a1*(a2prime-1) that the sparse transition matrix is built
%     from. A wrong stride there gives a valid but wrong index, so it would show up here as a
%     nonzero rather than as an error.
% The Policy comparisons matter as much as the V ones: with two endogenous states Policy carries a
% separate a2prime channel, so a mis-decoded second asset shows up there even if V looks fine.
n_a2A=[31,11]; % deliberately small, this is a correctness check and not a timing one
a1_grid2A=5*linspace(0,1,n_a2A(1))'.^3;
a2_grid2A=3*linspace(0,1,n_a2A(2))'.^3;
a_grid2A=[a1_grid2A; a2_grid2A]; % stacked, as the toolkit expects for multiple endogenous states
Params.r2=0.03; % return on the second asset (the first pays Params.r)
if prod(n_d)>0
    ReturnFn2A=@(d,a1prime,a2prime,a1,a2,z,r,r2,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz_with2A(d,a1prime,a2prime,a1,a2,z,r,r2,w,sigma,eta,varphi);
else
    ReturnFn2A=@(a1prime,a2prime,a1,a2,z,r,r2,w,sigma) ReturnFn_nod_z_noe_nosemiz_with2A(a1prime,a2prime,a1,a2,z,r,r2,w,sigma);
end

% Reference: postGI2A, greedy0/sparse0
vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.preGI=0;
[VG2A,PolicyG2A]=ValueFnIter_InfHorz(n_d,n_a2A,n_z,d_grid,a_grid2A,z_grid,pi_z,ReturnFn2A,Params,DF,[],vfo);

% preGI2A vs postGI2A
vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.preGI=1;
[Vpre2A,Policypre2A]=ValueFnIter_InfHorz(n_d,n_a2A,n_z,d_grid,a_grid2A,z_grid,pi_z,ReturnFn2A,Params,DF,[],vfo);
fprintf('GI 2A: preGI vs postGI,                V   should be ~0: %.3e \n',max(abs(VG2A(:)-Vpre2A(:))));
fprintf('GI 2A: preGI vs postGI,                Pol should be  0: %.3e \n',max(abs(PolicyG2A(:)-Policypre2A(:))));

% postGI2A, howardssparse=1, against the same reference
vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=1; vfo.preGI=0;
[Vsp2A,Policysp2A]=ValueFnIter_InfHorz(n_d,n_a2A,n_z,d_grid,a_grid2A,z_grid,pi_z,ReturnFn2A,Params,DF,[],vfo);
fprintf('GI 2A: sparse1 vs sparse0 (postGI),    V   should be ~0: %.3e \n',max(abs(VG2A(:)-Vsp2A(:))));
fprintf('GI 2A: sparse1 vs sparse0 (postGI),    Pol should be  0: %.3e \n',max(abs(PolicyG2A(:)-Policysp2A(:))));

% endogenousvfi=1 with two endogenous states (the GI2A endogenous twins): the span is over (a1,a2) within each z
vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.preGI=0; vfo.endogenousvfi=1;
[VE2A,PolicyE2A]=ValueFnIter_InfHorz(n_d,n_a2A,n_z,d_grid,a_grid2A,z_grid,pi_z,ReturnFn2A,Params,DF,[],vfo);
fprintf('GI 2A: endoVFI postGI vs ref,          V   should be ~0: %.3e \n',max(abs(VG2A(:)-VE2A(:))));
fprintf('GI 2A: endoVFI postGI vs ref,          Pol should be  0: %.3e \n',max(abs(PolicyG2A(:)-PolicyE2A(:))));
vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=1; vfo.preGI=0; vfo.endogenousvfi=1;
[VE2A,PolicyE2A]=ValueFnIter_InfHorz(n_d,n_a2A,n_z,d_grid,a_grid2A,z_grid,pi_z,ReturnFn2A,Params,DF,[],vfo);
fprintf('GI 2A: endoVFI sparse1 (postGI),       V   should be ~0: %.3e \n',max(abs(VG2A(:)-VE2A(:))));
fprintf('GI 2A: endoVFI sparse1 (postGI),       Pol should be  0: %.3e \n',max(abs(PolicyG2A(:)-PolicyE2A(:))));
vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.preGI=1; vfo.endogenousvfi=1;
[VE2A,PolicyE2A]=ValueFnIter_InfHorz(n_d,n_a2A,n_z,d_grid,a_grid2A,z_grid,pi_z,ReturnFn2A,Params,DF,[],vfo);
fprintf('GI 2A: endoVFI preGI vs postGI ref,    V   should be ~0: %.3e \n',max(abs(VG2A(:)-VE2A(:))));
fprintf('GI 2A: endoVFI preGI vs postGI ref,    Pol should be  0: %.3e \n',max(abs(PolicyG2A(:)-PolicyE2A(:))));
clear VE2A PolicyE2A

% preGI2A with howardssparse=1 -- NOT YET IMPLEMENTED (postGI only). Included, commented out:
% vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=1; vfo.preGI=1;
% [Vpre2Asp,Policypre2Asp]=ValueFnIter_InfHorz(n_d,n_a2A,n_z,d_grid,a_grid2A,z_grid,pi_z,ReturnFn2A,Params,DF,[],vfo);
% fprintf('GI 2A: preGI sparse1 vs postGI,        V   should be ~0: %.3e \n',max(abs(VG2A(:)-Vpre2Asp(:))));

clear VG PolicyG VpreA VpreB VpreC PolicypreA PolicypreB PolicypreC
clear VG2A PolicyG2A Vpre2A Policypre2A Vsp2A Policysp2A

%% Part 4: big n_a (=1500), with and without grid interpolation should give very similar V
% (Not exactly equal: GI refines the choice between grid points, so V_GI is a
%  touch higher; with 1500 points the difference is tiny.)
vfo=vfoptions; % without GI
[Vbig,~]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
vfo=vfoptions; vfo.gridinterplayer=1; vfo.ngridinterp=5; % with GI
[VbigGI,~]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('big n_a=%d: with vs without GI, V should be very similar (small): %.3e \n',n_a_big,max(abs(Vbig(:)-VbigGI(:))));

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
fprintf('postGI default window: |preGI-postGI| = %.3e,  min(Vpre-Vpost) = %.3e  (>=0 => preGI weakly higher, i.e. postGI under-searches) \n',max(abs(Vpre(:)-Vpost(:))),min(Vpre(:)-Vpost(:)));
mad=floor((n_a-1)/2); % widest VALID window: postGI clamps aprimeshifter to [1+mad, N_a-mad], so this covers the whole grid = full search
vfo=vfoGI; vfo.preGI=0; vfo.maxaprimediff=mad; [VpostW,~] =ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
fprintf('postGI WIDE window (maxaprimediff=%d): |preGI-postGI_wide| = %.3e  (if ~0 => the window was the cause; preGI is the exact reference) \n',mad,max(abs(Vpre(:)-VpostW(:))));

output=struct();

end
