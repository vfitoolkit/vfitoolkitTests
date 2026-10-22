function output=CoreInfHorzVFIAlgo_InfNaNcompare(label,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions,figure_c)
% Engine for the deliberate-trigger tests of the InfHorz grid interpolation layer's Howards step on a
% model where V contains -Inf (InfHorzGI_InfNaN_proposal.md at the toolkit root). Same role and
% same inputs as CoreInfHorzVFIAlgo_algocompare, for one model. Handles the plain models, semiz
% (vfoptions.n_semiz/semiz_grid/SemiExoStateFn, d2 the last d) and sophisticated quasi-hyperbolic
% (vfoptions.exoticpreferences) models.
%
% THE DEFECT. Every Howards step in the postGI raws evaluates the chosen fine-grid point as
%    lowerprob.*V(lowerind,:)+(1-lowerprob).*V(lowerind+1,:)
% then multiplies by pi_z and sets NaN to 0 (which is meant for pi_z==0 only). lowerprob is exactly
% 0 whenever the chosen point is a coarse node (the raws code that as lower=node-1 with all the
% weight on the upper node) and exactly 1 at the window's bottom node. A zero weight on a node where
% V is -Inf gives 0*(-Inf)=NaN, which wipes the whole (a,z,z') entry, and the isnan line then prices
% that continuation at 0 instead of V(node,z'). With u<0 that inflates V at every state next to a
% dead-end state. The fix, already in the _e raws, zeroes each weighted term before the sum. (The
% semiz raws combine EVz(EVindex) and EVz(EVindex+1) the same way before .*Pdz; the sophisticated
% quasi-hyperbolic raws combine the pi-summed EVlowerH/EVupperH with (1-wupperH)/wupperH.)
%
% WHY 'GI vs pure VFI' AT THE DEFAULTS CANNOT SEE IT. Howards is switched off after
% vfoptions.maxhowards (default 500) iterations of the fine-grid loop, and the plain Bellman
% iterations that follow converge from the inflated V to the correct fixed point. So at the defaults
% the converged V and Policy are right even on the defective raws; the only trace is a few hundred
% extra iterations. Every config is therefore solved twice:
%   (a) at the defaults: census of V, firing count, V against ValueFnFromPolicy on its own policy,
%       V against the pure-VFI GI solve. Expected to pass before AND after the fix; it is here so
%       that the laundering is on the record, and so the firing count proves the trigger fired.
%   (b) with Howards never switched off (maxhowards=maxiter=2000): on the defective raws the
%       alternation of a correct Bellman step and a wrong Howards step cannot converge, so the solve
%       hits maxiter (warning) and V is wrong; on the fixed raws it converges and matches pure VFI.
%       This is the red/green discriminator. (The endogenous twins can converge even on the
%       defective raws when the mispricing is the same for every state in an exogenous row, since
%       their span criterion ignores a within-row constant; the below-NBL models catch them.)
%
% THE FIRING COUNT is computed from the policy the solve returns, which is the (lower node, L2)
% pair the Howards step works from: L2==n2short+2 is a chosen coarse node coded as 'all weight on
% the upper node', so the zero weight is on the LOWER node; L2==1 is the window's bottom node, zero
% weight on the UPPER node. A state fires when its zero-weight node is -Inf at some (z',e') it can
% reach; with semiz, reachable under the d2 it chose. A model where this is 0 is vacuous for this
% defect, and the subcode says so.
%
% The greedy Howards configs (howardsgreedy=1,2,3) solve a linear system with -Inf in the return
% vector; the toolkit itself warns that greedy Howards cannot be used when V has -Inf. They are run
% so that the diary records what they do, but nothing they print is a check.
%
% Sophisticated quasi-hyperbolic is run at beta0=1 ONLY: at beta0<1 the continuation being
% interpolated is kinked and the bank documents that convergence itself is then a diagnostic, not
% a check (see QHDInfHorz_nod_z_GI), which would confound the red/green reading here. At beta0=1 the
% two quasi-hyperbolic postGI raws run the same Howards code on Vunderbar and the solution is exact.

fprintf('\n================ %s ================\n',label);
DF=DiscountFactorParamNames;

N_d=prod(n_d);
N_z=prod(n_z);
l_a=length(n_a);
n_a1=n_a(1);
N_arest=prod(n_a)/n_a1; % 1 unless two endogenous states
withe=isfield(vfoptions,'n_e');
withsemiz=isfield(vfoptions,'n_semiz');
withQH=isfield(vfoptions,'exoticpreferences');
if withe
    N_e=prod(vfoptions.n_e);
    pi_e_cpu=gather(vfoptions.pi_e(:)'); % [1,N_e]
else
    N_e=1;
    pi_e_cpu=1;
end
if N_z>0
    pi_z_cpu=gather(pi_z);
else
    pi_z_cpu=1;
end
N_zonly=max(N_z,1);
if withsemiz
    % The exogenous state is the joint (semiz,z), semiz fastest, and its transition depends on the chosen d2 (the last d).
    % pi_semiz is evaluated here on the CPU from the bank's SemiExoStateFn (its parameters are the bank setup's).
    N_semiz=prod(vfoptions.n_semiz);
    semiz_grid_cpu=gather(vfoptions.semiz_grid);
    N_d2=n_d(end);
    d2_grid_cpu=gather(d_grid(end-N_d2+1:end));
    pi_semiz_cpu=zeros(N_semiz,N_semiz,N_d2);
    for ss=1:N_semiz
        for sp=1:N_semiz
            for dd=1:N_d2
                pi_semiz_cpu(ss,sp,dd)=vfoptions.SemiExoStateFn(semiz_grid_cpu(ss),semiz_grid_cpu(sp),d2_grid_cpu(dd),Params.probfindjob,Params.problosejob);
            end
        end
    end
else
    N_semiz=1;
    N_d2=1;
end
N_zeff=N_semiz*N_zonly;
sz4=[n_a1,N_arest,N_zeff,N_e]; % V as [a1, rest of a, (semiz,)z, e], with singleton dims where the model has none
sz5=[sz4,N_d2]; % deadreach also by the chosen d2 (semiz only)
l_d=length(n_d)*(N_d>0);
n2short=5; % vfoptions.ngridinterp

%% No-GI reference: defines the dead-end set (the states where V is -Inf)
vfo=vfoptions; vfo.howardsgreedy=0; vfo.howardssparse=0;
[Vref,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
Vref=gather(Vref);
fprintf('noGI: V census: finite %i, -Inf %i, +Inf %i, NaN %i \n',nnz(isfinite(Vref)),nnz(Vref==-Inf),nnz(Vref==Inf),nnz(isnan(Vref)));
fprintf('noGI: NaN count,                                         should be  0: %i \n',nnz(isnan(Vref)));
fprintf('noGI: V has -Inf states (else this model is vacuous),    should be one: %i \n',nnz(Vref==-Inf)>0);
nminfref=nnz(Vref==-Inf);

Vref4=reshape(Vref,sz4);
Vd=(Vref4==-Inf);
finref=isfinite(Vref4);
% The proposal's policy-free stand-in: finite states with a -Inf neighbour on the a1 grid at the same (z,e)
lowerdead=false(sz4); lowerdead(2:end,:,:,:)=Vd(1:end-1,:,:,:);
upperdead=false(sz4); upperdead(1:end-1,:,:,:)=Vd(2:end,:,:,:);
fprintf('noGI: finite states with a -Inf lower neighbour (same z,e): %i, with a -Inf upper neighbour: %i \n',nnz(finref&lowerdead),nnz(finref&upperdead));
fprintf('noGI: some finite state has a -Inf neighbour,            should be one: %i \n',(nnz(finref&lowerdead)+nnz(finref&upperdead))>0);

% The plain (no-GI) sparse Howards raws multiply a sparse matrix into a vector that holds -Inf too, and GPU sparse*vector gives NaN
% on a -Inf (the probe at the top of the Part), so they are checked here against the indexed reference (class B of the proposal).
% (Not for quasi-hyperbolic, which has no sparse raws.)
if ~withQH
    vfo=vfoptions; vfo.howardsgreedy=0; vfo.howardssparse=1;
    [Vrefsp,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    Vrefsp=gather(Vrefsp);
    fprintf('noGI sparse1: V census: finite %i, -Inf %i, +Inf %i, NaN %i \n',nnz(isfinite(Vrefsp)),nnz(Vrefsp==-Inf),nnz(Vrefsp==Inf),nnz(isnan(Vrefsp)));
    fprintf('noGI sparse1: -Inf count minus noGI -Inf count,              should be  0: %i \n',nnz(Vrefsp==-Inf)-nminfref);
    finboth=isfinite(Vrefsp)&isfinite(Vref);
    fprintf('noGI sparse1 vs noGI sparse0 on the finite set,         V   should be ~0: %.3e \n',max(abs(Vrefsp(finboth)-Vref(finboth))));
    vfo.maxiter=2000; vfo.maxhowards=2000;
    lastwarn('');
    [Vrefsp,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    msgsp=lastwarn;
    Vrefsp=gather(Vrefsp);
    fprintf('noGI sparse1, Howards never off: converged (no maxiter warning), should be one: %i \n',~contains(msgsp,'maximum number of iterations'));
    finboth=isfinite(Vrefsp)&isfinite(Vref);
    fprintf('noGI sparse1, Howards never off: vs noGI sparse0 on the finite set, V should be ~0: %.3e \n',max(abs(Vrefsp(finboth)-Vref(finboth))));
end

% deadreach(a1,arest,zeff,e,d2): is V(a1,arest,zeff',e') -Inf at some (zeff',e') reachable from (zeff,e) under d2?
% (a Howards weight-0 term on that node, from a state at (zeff,e) that chose d2, is what fires the defect)
deadreach=false(sz5);
Vd2=reshape(Vd,[n_a1*N_arest,N_zeff*N_e]);
for zz=1:N_zeff
    for ee=1:N_e
        for dd=1:N_d2
            if withsemiz
                ss=rem(zz-1,N_semiz)+1; % semiz index of the joint state
                zo=ceil(zz/N_semiz); % z index of the joint state
                reachz=(pi_semiz_cpu(ss,:,dd)'>0)&(pi_z_cpu(zo,:)>0); % [N_semiz,N_zonly]: joint (semiz',z') reachable, semiz' fastest
                reach=(reachz(:))&(pi_e_cpu>0); % [N_zeff,N_e]
            else
                reach=(pi_z_cpu(zz,:)'>0)&(pi_e_cpu>0); % [N_zeff,N_e], iid e so the e row is the same for every e
            end
            deadreach(:,:,zz,ee,dd)=reshape(any(Vd2(:,reach(:)),2),[n_a1,N_arest]);
        end
    end
end

%% Pure VFI with the grid interpolation layer: no Howards on either grid, so no defect, the reference for the GI solves
vfoGI=vfoptions; vfoGI.gridinterplayer=1; vfoGI.ngridinterp=n2short;
if withQH
    % Sophisticated quasi-hyperbolic with the grid interpolation layer has no default window (the dispatcher errors); the
    % value its error message recommends. And no sparse Howards: the quasi-hyperbolic postGI raws have no sparse variant.
    vfoGI.maxaprimediff=ceil(n_a1/(Params.beta0*4));
    vfoGI.howardssparse=0;
end
vfo=vfoGI; vfo.howardsgreedy=0; vfo.howardssparse=0; vfo.howards=0;
[Vpure,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
Vpure=gather(Vpure);
fprintf('GI pure VFI (howards=0): V census: finite %i, -Inf %i, +Inf %i, NaN %i \n',nnz(isfinite(Vpure)),nnz(Vpure==-Inf),nnz(Vpure==Inf),nnz(isnan(Vpure)));
fprintf('GI pure VFI (howards=0): NaN count,                      should be  0: %i \n',nnz(isnan(Vpure)));
fprintf('GI pure VFI (howards=0): -Inf count minus noGI -Inf count, should be  0: %i \n',nnz(Vpure==-Inf)-nminfref);
% (GI and no-GI V differ genuinely, interpolation refines the choice, so no identity between them is printed; the figure shows both)
finboth=isfinite(Vpure)&isfinite(Vref);
fprintf('GI pure VFI vs noGI on the joint finite set, max abs diff (diagnostic, interpolation refines the choice): %.3e \n',max(abs(Vpure(finboth)-Vref(finboth))));

%% The configs. {name, howardsgreedy, howardssparse, postGIrepeat, kind, endogenousvfi}
% Without z and without e the grid interpolation layer is a single postGI solve that reads neither howardsgreedy nor
% howardssparse (see algocompare), so those models get the reference solve and the postGIrepeat layer only.
% With two endogenous states howardsgreedy>0 is refused by the toolkit. With e, HowardMix (2,3) is refused.
% With semiz the grid interpolation layer is postGI only: no postGIrepeat, and sparse and greedy Howards only without a d1
% (with a d1 howardssparse is ignored, so comparing it would compare a solve with itself).
% Quasi-hyperbolic: one postGI raw per d-ness, no sparse, no greedy, no endogenous VFI (refused for exotic preferences).
configs=cell(0,5);
configs(end+1,:)={'greedy0/sparse0        ',0,0,0,'check'};
if ~withsemiz && ~withQH
    configs(end+1,:)={'greedy0/sparse0/repeat1',0,0,1,'check'};
end
if withsemiz
    if length(n_d)==1
        configs(end+1,:)={'greedy0/sparse1        ',0,1,0,'check'};
        configs(end+1,:)={'greedy1 (HowardGreedy) ',1,0,0,'info'};
    end
elseif ~withQH
    if N_z>0 || withe
        configs(end+1,:)={'greedy0/sparse1        ',0,1,0,'check'};
        configs(end+1,:)={'greedy0/sparse1/repeat1',0,1,1,'check'};
    end
    if N_z>0 && l_a==1
        configs(end+1,:)={'greedy1 (HowardGreedy) ',1,0,0,'info'};
        if ~withe
            configs(end+1,:)={'greedy2 (HowardMix)    ',2,0,0,'info'};
            configs(end+1,:)={'greedy3 (HowardMix2)   ',3,0,0,'info'};
        end
    end
end
% Every check config gets its vfoptions.endogenousvfi=1 twin: the EndogenousVFI raws (RelativeVFI without z, and the semiz
% twins of both) have their own copies of the Howards step, so they fire the same defect; the e twins already carry the
% per-term form (a positive control, like the e raws above).
nconfigs=size(configs,1);
configs(:,6)={0};
if ~withQH
    for cc=1:nconfigs
        if strcmp(configs{cc,5},'check')
            configs(end+1,:)={['endoVFI ',configs{cc,1}],configs{cc,2},configs{cc,3},configs{cc,4},'check',1};
        end
    end
end

for cc=1:size(configs,1)
    cname=configs{cc,1};
    vfo=vfoGI; vfo.howardsgreedy=configs{cc,2}; vfo.howardssparse=configs{cc,3}; vfo.postGIrepeat=configs{cc,4}; vfo.endogenousvfi=configs{cc,6};

    if strcmp(configs{cc,5},'info')
        % Greedy Howards: not a check, see the header. try/catch so that an error here does not kill the part.
        lastwarn('');
        try
            [VG,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
            VG=gather(VG);
            finboth=isfinite(VG)&isfinite(Vpure);
            fprintf('GI %s (informational, greedy Howards is not supported with -Inf in V): V census: finite %i, -Inf %i, +Inf %i, NaN %i; max abs diff vs pure VFI on the joint finite set %.3e; last warning: %s \n',...
                cname,nnz(isfinite(VG)),nnz(VG==-Inf),nnz(VG==Inf),nnz(isnan(VG)),max(abs(VG(finboth)-Vpure(finboth))),lastwarn);
        catch ME
            fprintf('GI %s (informational, greedy Howards is not supported with -Inf in V): errored: %s \n',cname,ME.message);
        end
        continue
    end

    %% (a) at the defaults
    [VG,PolicyG]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    VGcpu=gather(VG);
    fprintf('GI %s: V census: finite %i, -Inf %i, +Inf %i, NaN %i \n',cname,nnz(isfinite(VGcpu)),nnz(VGcpu==-Inf),nnz(VGcpu==Inf),nnz(isnan(VGcpu)));
    fprintf('GI %s: NaN count,                                 should be  0: %i \n',cname,nnz(isnan(VGcpu)));
    fprintf('GI %s: -Inf count minus noGI -Inf count,          should be  0: %i \n',cname,nnz(VGcpu==-Inf)-nminfref);
    fprintf('GI %s: finite set differs from pure VFI (count),  should be  0: %i \n',cname,nnz(isfinite(VGcpu)~=isfinite(Vpure)));
    finboth=isfinite(VGcpu)&isfinite(Vpure);
    fprintf('GI %s: vs pure VFI on the finite set,        V   should be ~0: %.3e \n',cname,max(abs(VGcpu(finboth)-Vpure(finboth))));

    % The firing count, from the returned policy (see the header)
    PolicyG2=reshape(gather(PolicyG),[size(PolicyG,1),prod(sz4)]);
    L1=PolicyG2(l_d+1,:); % lower node on the a1 grid
    L2=PolicyG2(end-1,:); % 1 (all weight on L1) to n2short+2 (all weight on L1+1); the last row is the L2flag
    if l_a>=2
        iarp=PolicyG2(l_d+2,:); % next period's index on the rest of the endogenous states
    else
        iarp=ones(1,prod(sz4));
    end
    if withsemiz
        id2=PolicyG2(l_d,:); % the chosen d2 (the last d) sets which semiz' are reachable
    else
        id2=ones(1,prod(sz4));
    end
    [~,~,iz,ie]=ind2sub(sz4,1:prod(sz4));
    drlower=reshape(deadreach(sub2ind(sz5,L1,iarp,iz,ie,id2)),[1,prod(sz4)]); % reshape: when deadreach is a column (no z, no e, one endogenous state) indexing it with a row gives a column
    drupper=reshape(deadreach(sub2ind(sz5,L1+1,iarp,iz,ie,id2)),[1,prod(sz4)]);
    finG=isfinite(VGcpu(:))';
    fireslower=nnz(finG & L2==n2short+2 & drlower);
    firesupper=nnz(finG & L2==1 & drupper);
    fprintf('GI %s: Howards puts a zero weight on a -Inf node at %i finite states (dead lower node %i, dead upper node %i) \n',cname,fireslower+firesupper,fireslower,firesupper);
    fprintf('GI %s: the trigger fired,                         should be one: %i \n',cname,(fireslower+firesupper)>0);

    % ValueFnFromPolicy on its own policy (it already carries the per-term form)
    Vfp=gather(ValueFnFromPolicy_InfHorz(PolicyG,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfo));
    fprintf('GI %s: finite set differs from ValueFnFromPolicy (count), should be  0: %i \n',cname,nnz(isfinite(VGcpu)~=isfinite(Vfp)));
    finboth=isfinite(VGcpu)&isfinite(Vfp);
    fprintf('GI %s: vs ValueFnFromPolicy on the finite set, V should be ~0: %.3e \n',cname,max(abs(VGcpu(finboth)-Vfp(finboth))));

    %% (b) Howards never switched off
    vfo.maxiter=2000; vfo.maxhowards=2000;
    lastwarn('');
    [VGh,PolicyGh]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    msgh=lastwarn;
    VGhcpu=gather(VGh);
    fprintf('GI %s, Howards never off: converged (no maxiter warning), should be one: %i \n',cname,~contains(msgh,'maximum number of iterations'));
    fprintf('GI %s, Howards never off: V census: finite %i, -Inf %i, +Inf %i, NaN %i \n',cname,nnz(isfinite(VGhcpu)),nnz(VGhcpu==-Inf),nnz(VGhcpu==Inf),nnz(isnan(VGhcpu)));
    fprintf('GI %s, Howards never off: NaN count,                        should be  0: %i \n',cname,nnz(isnan(VGhcpu)));
    fprintf('GI %s, Howards never off: finite set differs from pure VFI (count), should be  0: %i \n',cname,nnz(isfinite(VGhcpu)~=isfinite(Vpure)));
    finboth=isfinite(VGhcpu)&isfinite(Vpure);
    fprintf('GI %s, Howards never off: vs pure VFI on the finite set, V should be ~0: %.3e \n',cname,max(abs(VGhcpu(finboth)-Vpure(finboth))));
    Vfph=gather(ValueFnFromPolicy_InfHorz(PolicyGh,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfo));
    finboth=isfinite(VGhcpu)&isfinite(Vfph);
    fprintf('GI %s, Howards never off: vs ValueFnFromPolicy on the finite set, V should be ~0: %.3e \n',cname,max(abs(VGhcpu(finboth)-Vfph(finboth))));

    if cc==1 % keep the reference config for the figure and the output
        VGdef=VGcpu;
        VGnever=VGhcpu;
        output.fireslower=fireslower;
        output.firesupper=firesupper;
    end
end

%% Figure: V(a1) at the first and last exogenous state, over the first nodes of the a1 grid (where the dead region and the
% constrained states sit), no GI against the GI solves. -Inf does not plot, so the dead region shows as a gap.
Vpure4=reshape(Vpure,sz4);
VGdef4=reshape(VGdef,sz4);
VGnever4=reshape(VGnever,sz4);
a1_grid=gather(a_grid(1:n_a1));
na1plot=min(n_a1,30);
zplot=[1,N_zeff];
figure(figure_c)
for pp=1:2
    subplot(1,2,pp)
    plot(a1_grid(1:na1plot),Vref4(1:na1plot,1,zplot(pp),1),'k-o',...
        a1_grid(1:na1plot),Vpure4(1:na1plot,1,zplot(pp),1),'g-',...
        a1_grid(1:na1plot),VGdef4(1:na1plot,1,zplot(pp),1),'b--',...
        a1_grid(1:na1plot),VGnever4(1:na1plot,1,zplot(pp),1),'r:')
    legend('no GI','GI pure VFI','GI defaults','GI Howards never off','Location','southeast')
    title([label,': exogenous state index ',num2str(zplot(pp))])
    xlabel(['a (first ',num2str(na1plot),' nodes)'])
    ylabel('V')
end

end
