function output=CoreInfHorzVFIAlgo_ScanHowardsSettings(Params,DiscountFactorParamNames)
% SCAN 3: the Howards accelerators -- do they give the same answer as no Howards, and are they
% actually faster? Done for the z (no e) models, with and without d, at n_z = 5, 15, 25, 75 and
% n_a = 100, 200, 500. (n_z is swept because the vfoptions.howardssparse default triggers on
%  N_z>100, but the runtimes had only ever been measured at n_z=5.)
%
% NOTE: vfoptions.maxhowards=0 is the off-switch for Howards. vfoptions.howards=0 only turns off the
% ITERATED Howards; the greedy paths never use vfoptions.howards (they are gated by maxhowards alone),
% so howards=0 would be a no-op for any howardsgreedy>0 and the comparison would be vacuous.
%
% Why the timings matter: Howards is only an accelerator, so a wrong Howards step still converges to
% the correct V (the outer value function iteration fixes it up). Hence the V/Policy checks alone
% CANNOT detect a broken accelerator -- only the runtime can.
%
% For the ITERATED Howards we try howards=40,80,120 (rather than the default of 150), so that the
% iterated and the greedy Howards runtimes can be compared against each other. The greedy configs
% get a single run, as they ignore vfoptions.howards.
%
% The GI sparse config is run at lowmemory=0 and lowmemory=1. Until recently sparse+GI was reachable
% only at lowmemory=1, and the raw behind it did not actually loop over z, so the old lowmem1
% timings were really lowmemory=0 timings. Both are now measured, so the trade-off is visible.
%
% Speedups are printed rather than asserted; a value below 1 means that accelerator is not earning
% its keep (or is broken).
%
% Greedy Howards (howardsgreedy=1) is skipped whenever n_a=500 or n_z=75: it solves a sparse linear
% system of size N_a*N_z, and is already clearly slower than iterated Howards at grid sizes well
% below these, so those cells print 'not run as slow' in place of a runtime.

DF=DiscountFactorParamNames;

n_z_list=[5,15,25,75];
n_a_list=[100,200,500];
howardslist=[40,80,120];

% Each config is {label, howardsgreedy, howardssparse, lowmemory, sweep howards?}
cfg_noGI={{'greedy0/sparse0     ',0,0,0,1},...
          {'greedy0/sparse1     ',0,1,0,1},...
          {'lowmem1/sparse0     ',0,0,1,1},...
          {'lowmem1/sparse1     ',0,1,1,1},...
          {'greedy1 (greedy)    ',1,0,0,0}};
% Note: sparse1 is run at both lowmemory settings. With d they differ in how the two refined return
% matrices are built (one z at a time, vs whole), so lowmem1 should be slower but lighter on memory;
% without d lowmemory has no effect on this code path and the two should time the same.
cfg_GI  ={{'greedy0/sparse0     ',0,0,0,1},...
          {'sparse1/lowmem0     ',0,1,0,1},...
          {'sparse1/lowmem1     ',0,1,1,1},...
          {'greedy2 (HowardMix) ',2,0,0,1},...
          {'greedy3 (HowardMix2)',3,0,0,1},...
          {'greedy1 (greedy)    ',1,0,0,0}};

fprintf('\n================ SCAN 3: Howards accelerators (same answer, and faster?) ================\n');

for zz=1:length(n_z_list)
    n_z=n_z_list(zz);
    [z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.03,n_z); % same discretization as the bank setup
    z_grid=exp(z_grid);

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

        for aa=1:length(n_a_list)
            n_a=n_a_list(aa);
            a_grid=5*linspace(0,1,n_a)'.^3;
            fprintf('\n---- %s, n_z=%d, n_a=%d ----\n',modelstr,n_z,n_a);

            vfoptions=struct(); vfoptions.verbose_advice=0;
            vfoGI=vfoptions; vfoGI.gridinterplayer=1; vfoGI.ngridinterp=5;

            % Warm up (so the first timed solve is not paying GPU/JIT startup costs)
            [~,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions);

            for gi=0:1
                if gi==0
                    vfobase=vfoptions; cfg=cfg_noGI; gistr='noGI';
                else
                    vfobase=vfoGI; cfg=cfg_GI; gistr='GI  ';
                end

                % Howards OFF: pure value function iteration, no acceleration of any kind
                vfo=vfobase; vfo.maxhowards=0;
                tic; [Voff,Policyoff]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); wait(gpuDevice); tOFF=toc;
                fprintf('%s Howards OFF (maxhowards=0), pure VFI runtime: %2.4f seconds \n',gistr,tOFF);

                for cc=1:length(cfg)
                    if gi==0 && cfg{cc}{4}==1 && n_d==0
                    % Without GI, lowmemory=1 exists only for refinement, which needs a d. A model
                    % with no d has no lowmemory option at all, so this config does not exist here.
                    fprintf('%s %s             not run, no lowmemory without d \n',gistr,cfg{cc}{1});
                    continue
                end
                if cfg{cc}{2}==1 && (n_a==500 || n_z==75)
                        % Greedy Howards solves a sparse linear system of size N_a*N_z, and the smaller grid
                        % sizes already show it losing badly to iterated Howards, so there is nothing to learn
                        % from paying for it here.
                        fprintf('%s %s             vs OFF, speedup >1: not run as slow \n',gistr,cfg{cc}{1});
                        continue
                    end
                    if cfg{cc}{5}==1 % iterated Howards, so sweep the number of Howards iterations
                        for hh=1:length(howardslist)
                            vfo=vfobase; vfo.howardsgreedy=cfg{cc}{2}; vfo.howardssparse=cfg{cc}{3}; vfo.lowmemory=cfg{cc}{4}; vfo.howards=howardslist(hh);
                            tic; [Vc,Policyc]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); wait(gpuDevice); tc=toc;
                            fprintf('%s %s howards=%3d vs OFF, V ~0: %2.8f, Pol 0: %2.8f, speedup >1: %2.2f \n',gistr,cfg{cc}{1},howardslist(hh),max(abs(Voff(:)-Vc(:))),max(abs(Policyoff(:)-Policyc(:))),tOFF/tc);
                        end
                    else % greedy Howards, which ignores vfoptions.howards
                        vfo=vfobase; vfo.howardsgreedy=cfg{cc}{2}; vfo.howardssparse=cfg{cc}{3}; vfo.lowmemory=cfg{cc}{4};
                        tic; [Vc,Policyc]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); wait(gpuDevice); tc=toc;
                        fprintf('%s %s             vs OFF, V ~0: %2.8f, Pol 0: %2.8f, speedup >1: %2.2f \n',gistr,cfg{cc}{1},max(abs(Voff(:)-Vc(:))),max(abs(Policyoff(:)-Policyc(:))),tOFF/tc);
                    end
                end
            end
        end
    end

end % n_z loop

output=struct();

end
