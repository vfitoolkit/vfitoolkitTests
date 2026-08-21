function output=CoreInfHorzVFIAlgo_ScanHowardsSettingsWith2A(Params,DiscountFactorParamNames)
% SCAN 5: the Howards settings for models with TWO endogenous states.
%
% Three parts:
%   Part A: without the grid interpolation layer, sweep the howardssparse/lowmemory configs.
%           With two endogenous states the pure discretization raws just see N_a=prod(n_a), so this is
%           the same code as the one-endogenous-state scan at a larger N_a. It is included because two
%           endogenous states is how models actually reach a large N_a (fig 17 of CoreInfHorzTests is
%           101-by-51=5151), so this is the realistic setting for the N_a part of the
%           vfoptions.howardssparse default, and the only place where a big N_a comes with a small
%           n_z. That combination is what showed the default has to be set on N_z rather than on
%           N_a*N_z: two endogenous state cells with N_a over 10000 and n_z of 5 or 25 are ties.
%   Part B: with the grid interpolation layer, compare howardssparse=0 against howardssparse=1, and
%           sweep the number of Howards iterations. howardsgreedy>0 is a deliberate error for two
%           endogenous states so it is not swept. The sparse raws here are new, and the non-sparse
%           ones only recently got a working Howards in their fine-grid stage at all, so both the
%           agreement between them and their relative runtimes are worth seeing.
%   Part C: check that the defaults actually dispatch. Solve with nothing set but gridinterplayer and
%           ngridinterp, and confirm it both runs and matches the explicit reference. The defaults are
%           size dependent, and with two endogenous states N_a is large enough to trip size-based
%           rules, so this is a real check rather than a formality.
%
% These grids can be big, so every solve is wrapped in try/catch and a failure prints in place of the
% runtime rather than killing the script.

DF=DiscountFactorParamNames;

Params.r2=0.03; % return on the second asset (the first pays Params.r)

% (n_z, n_a) pairs, where n_a is now a pair. Not a full cross product.
% n_z=50 is in here because it is exactly the threshold in the vfoptions.howardssparse default
% ('N_z>=50 && N_a>=500' without the grid interpolation layer). Every other two endogenous state
% cell has n_z of 5 or 25, so that branch of the default has never been reached with two endogenous
% states, and it is the one place the right setting could plausibly differ from the one endogenous
% state case. These two cells are the biggest in the scan (N_a*N_z of 258k and 536k, against 225k
% for the largest cell of the big grids scan), so expect them to be slow and possibly to run out of
% memory for the d model.
cells={{5,[51,21]},{5,[101,51]},{25,[101,51]},{50,[101,51]},{5,[151,71]},{25,[151,71]},{50,[151,71]}};
howardslist=[40,80,120];

% Each config is {label, howardsgreedy, howardssparse, lowmemory, needs d?}
cfg_noGI={{'greedy0/sparse0     ',0,0,0,0},...
          {'greedy0/sparse1     ',0,1,0,0},...
          {'lowmem1/sparse0     ',0,0,1,1},...
          {'lowmem1/sparse1     ',0,1,1,1}};
% With two endogenous states the grid interpolation layer implements greedy0 with either
% howardssparse setting (howardsgreedy>0 is a deliberate error there, so it is not swept)
cfg_GI  ={{'greedy0/sparse0     ',0,0,0,0},...
          {'sparse1/lowmem0     ',0,1,0,0}};

fprintf('\n================ SCAN 5: Howards settings, two endogenous states ================\n');

for cc_c=1:length(cells)
n_z=cells{cc_c}{1};
n_a=cells{cc_c}{2};
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.03,n_z); % same discretization as the bank setup
z_grid=exp(z_grid);
a1_grid=5*linspace(0,1,n_a(1))'.^3;
a2_grid=3*linspace(0,1,n_a(2))'.^3;
a_grid=[a1_grid; a2_grid]; % stacked, as the toolkit expects for multiple endogenous states

for mm=1:2
    if mm==1
        n_d=0; d_grid=[];
        ReturnFn=@(a1prime,a2prime,a1,a2,z,r,r2,w,sigma) ReturnFn_nod_z_noe_nosemiz_with2A(a1prime,a2prime,a1,a2,z,r,r2,w,sigma);
        modelstr='nod, z, noe, 2A';
    else
        n_d=9; d_grid=linspace(0,1,n_d)';
        ReturnFn=@(d,a1prime,a2prime,a1,a2,z,r,r2,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz_with2A(d,a1prime,a2prime,a1,a2,z,r,r2,w,sigma,eta,varphi);
        modelstr='d, z, noe, 2A';
    end

    fprintf('\n---- %s, n_z=%d, n_a=[%d,%d] (N_a=%d, N_a*N_z=%d) ----\n',...
        modelstr,n_z,n_a(1),n_a(2),prod(n_a),prod(n_a)*n_z);

    vfoptions=struct(); vfoptions.verbose_advice=0;
    vfoGI=vfoptions; vfoGI.gridinterplayer=1; vfoGI.ngridinterp=5;

    for gi=0:1
        if gi==0
            vfobase=vfoptions; cfg=cfg_noGI; gistr='noGI'; % Part A
        else
            vfobase=vfoGI; cfg=cfg_GI; gistr='GI  '; % Part B
        end

        % Howards OFF reference: pure value function iteration, no acceleration of any kind
        Voff=[]; Policyoff=[]; tOFF=NaN;
        vfo=vfobase; vfo.maxhowards=0;
        try
            tic; [Voff,Policyoff]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); wait(gpuDevice); tOFF=toc;
            fprintf('%s Howards OFF (maxhowards=0), pure VFI runtime: %2.4f seconds \n',gistr,tOFF);
            % Pull the reference back to the cpu, so a later failure cannot leave it stale.
            % V and Policy are small, so this costs nothing.
            Voff=gather(Voff); Policyoff=gather(Policyoff);
        catch ME
            fprintf('%s Howards OFF (maxhowards=0), not run: %s \n',gistr,ME.message(1:min(end,90)));
        end

        for cc=1:length(cfg)
            if cfg{cc}{5}==1 && n_d==0
                % Without the grid interpolation layer, lowmemory=1 exists only for refinement, which needs a d
                fprintf('%s %s             not run, no lowmemory without d \n',gistr,cfg{cc}{1});
                continue
            end
            for hh=1:length(howardslist)
                vfo=vfobase; vfo.howardsgreedy=cfg{cc}{2}; vfo.howardssparse=cfg{cc}{3}; vfo.lowmemory=cfg{cc}{4}; vfo.howards=howardslist(hh);
                try
                    tic; [Vc,Policyc]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); wait(gpuDevice); tc=toc;
                    if isempty(Voff)
                        fprintf('%s %s howards=%3d runtime: %2.4f seconds (no OFF reference to compare against) \n',gistr,cfg{cc}{1},howardslist(hh),tc);
                    else
                        fprintf('%s %s howards=%3d vs OFF, V ~0: %2.8f, Pol 0: %2.8f, speedup >1: %2.2f \n',gistr,cfg{cc}{1},howardslist(hh),max(abs(Voff(:)-Vc(:))),max(abs(Policyoff(:)-Policyc(:))),tOFF/tc);
                    end
                catch ME
                    fprintf('%s %s howards=%3d not run: %s \n',gistr,cfg{cc}{1},howardslist(hh),ME.message(1:min(end,90)));
                end
            end
        end

        %% Part C: do the defaults dispatch, and do they give the reference answer?
        % Nothing is set here except the grid interpolation layer itself, so every Howards option
        % takes whatever ValueFnIter_InfHorz picks. With two endogenous states the grid interpolation
        % layer only implements howardssparse=0, so a default that resolved to 1 would error here.
        vfo=vfobase; % gridinterplayer/ngridinterp only, no Howards options set at all
        try
            tic; [Vd,Policyd]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); wait(gpuDevice); td=toc;
            if isempty(Voff)
                fprintf('%s defaults (nothing set)   runtime: %2.4f seconds (no OFF reference to compare against) \n',gistr,td);
            else
                fprintf('%s defaults (nothing set)   vs OFF, V ~0: %2.8f, Pol 0: %2.8f, speedup >1: %2.2f \n',gistr,max(abs(Voff(:)-Vd(:))),max(abs(Policyoff(:)-Policyd(:))),tOFF/td);
            end
        catch ME
            fprintf('%s defaults (nothing set)   NOT RUN, the defaults do not dispatch here: %s \n',gistr,ME.message(1:min(end,90)));
        end
    end
end

end % cells loop

output=struct();

end
