function output=CoreInfHorzVFIAlgo_ScanHowardsBigGrids(Params,DiscountFactorParamNames)
% SCAN 4: the Howards defaults on big grids.
%
% Purpose is narrow: decide whether the vfoptions.howardssparse default is set on the right rule.
% Without the grid interpolation layer that default is currently
%     if N_a>1200 && N_z>100, howardssparse=1, else howardssparse=0
% and CoreInfHorzVFIAlgo_ScanHowardsSettings never once fired it (it stops at n_a=500, n_z=75), so
% the rule has never actually been tested. The five (n_z,n_a) cells here straddle its two conditions:
%
%   (n_z,n_a)      N_a>1200   N_z>100    what it tells us
%   (75,1000)      no         no         neighbour, for the trend
%   (150,500)      no         yes        is N_z on its own enough?
%   (150,1000)     no         yes        interior
%   (75,1500)      yes        no         is N_a on its own enough?
%   (150,1500)     yes        yes        the only cell where the current rule fires
%
% Deliberately much smaller than SCAN 3, because most of what SCAN 3 checks is already settled and
% just slow to redo at these sizes:
%   - howardsgreedy=1,2,3 are not run at all. They were 0.01 to 0.25 by n_z=25 and get worse with
%     size, so vfoptions.howardsgreedy=0 is not in question.
%   - howards=120 is not run. It won 2 of the 120 cells in SCAN 3.
%
% These grids are big enough to run out of GPU memory, so every solve is wrapped in try/catch and a
% failure prints in place of the runtime rather than killing the script. In the first run of this
% scan every single d model GI cell ran out of memory, at both of the GI configs, so each config is
% now run at lowmemory=0 and lowmemory=1. The lowmemory=0 ones are expected to keep failing until
% there is a bigger GPU; the point is that lowmemory=1 loops over z when building the refined return
% matrices and so may fit, which would at least show what these grids can do.

DF=DiscountFactorParamNames;

% (n_z,n_a) pairs, not a full cross product
cells={[75,1000],[150,500],[150,1000],[75,1500],[150,1500]};
howardslist=[40,80];

% Each config is {label, howardsgreedy, howardssparse, lowmemory}
cfg_noGI={{'greedy0/sparse0     ',0,0,0},...
          {'greedy0/sparse1     ',0,1,0},...
          {'lowmem1/sparse0     ',0,0,1},...
          {'lowmem1/sparse1     ',0,1,1}};
% Both GI configs are run at each lowmemory setting. With a d, lowmemory=1 makes refinement build
% the two refined return matrices one z at a time rather than whole, which is the only way these
% cells fit: at lowmemory=0 the [N_d,N_aprime,N_a,N_z] array before the max over d is tens of GB.
cfg_GI  ={{'greedy0/sparse0/lm0 ',0,0,0},...
          {'greedy0/sparse0/lm1 ',0,0,1},...
          {'sparse1/lowmem0     ',0,1,0},...
          {'sparse1/lowmem1     ',0,1,1}};

fprintf('\n================ SCAN 4: Howards defaults on big grids ================\n');

for cc_c=1:length(cells)
n_z=cells{cc_c}(1);
n_a=cells{cc_c}(2);
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.03,n_z); % same discretization as the bank setup
z_grid=exp(z_grid);
a_grid=5*linspace(0,1,n_a)'.^3;

for mm=1:2
    if mm==1
        n_d=0; d_grid=[];
        ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
        modelstr='nod, z, noe';
    else
        n_d=9; d_grid=linspace(0,1,n_d)';
        ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
        modelstr='d, z, noe';
    end

    fprintf('\n---- %s, n_z=%d, n_a=%d (N_a*N_z=%d; current rule would set howardssparse=%d) ----\n',...
        modelstr,n_z,n_a,n_a*n_z,(n_a>1200 && n_z>100));

    vfoptions=struct(); vfoptions.verbose_advice=0;
    vfoGI=vfoptions; vfoGI.gridinterplayer=1; vfoGI.ngridinterp=5;

    for gi=0:1
        if gi==0
            vfobase=vfoptions; cfg=cfg_noGI; gistr='noGI';
        else
            vfobase=vfoGI; cfg=cfg_GI; gistr='GI  ';
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
            if mm==1 && cfg{cc}{4}==1
                % lowmemory=1 only ever means 'refinement builds the return matrix one z at a time',
                % and refinement needs a d. Without one there is nothing for it to do: without GI it
                % is an error, with GI it is silently ignored (SCAN 3 timed the two the same, 0.94 to
                % 1.01, across all twelve nod cells), so there is no point running these twice.
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
    end
end

end % cells loop

output=struct();

end
