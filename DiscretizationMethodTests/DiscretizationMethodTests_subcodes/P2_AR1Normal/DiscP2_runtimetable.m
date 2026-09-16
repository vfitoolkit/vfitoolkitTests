function DiscP2_runtimetable(outputP2,calib,znums)
% Runtime table for P2 (AR(1) with normal innovations)
%
% Reported, never asserted. Measured on the same calls the accuracy sweep read, so the runtime
% table and the accuracy figures share an x-axis by construction.
%
% This is the block where the runtime column starts mattering for method choice: Farmer-Toda
% solves an entropy problem per grid point, so it is orders of magnitude slower than the
% closed-form methods, and P2's accuracy figures are what say whether that buys anything.

fprintf('\n========== Runtimes for discretizing an AR(1) with normal innovations ========== \n')
v=version;
try
    g=gpuDevice; gname=g.Name;
catch
    gname='no gpu';
end
fprintf('(median of up to %i timed calls after 1 discarded warm-up; MATLAB %s; GPU: %s) \n',calib.nreps,v,gname)

nmlist={'discretizeAR1_Tauchen','discretizeAR1_Rouwenhorst','discretizeAR1_TauchenHussey','discretizeAR1_FarmerToda'};
fldlist={'Tauchen','Rouwenhorst','TauchenHussey','FarmerToda'};

for cal_c=1:3
    cname=calib.names{cal_c};
    fprintf('\n--- calibration %s (rho=%g) --- \n',cname,calib.(cname).rho);
    fprintf('%-34s',' ');
    for c_c=1:length(znums)
        fprintf('%11s',['znum=',num2str(znums(c_c))]);
    end
    fprintf(' \n');
    for m_c=1:4
        fprintf('%-34s',nmlist{m_c});
        t=outputP2.(fldlist{m_c}).(cname).runtime;
        for c_c=1:length(znums)
            fprintf('%11.6f',t(c_c));
        end
        fprintf(' \n');
    end
end

% The comparison worth reading off the table
tR=outputP2.Rouwenhorst.moderate.runtime(end);
tF=outputP2.FarmerToda.moderate.runtime(end);
fprintf('\nAt znum=%i on the moderate calibration, FarmerToda takes %2.1fx as long as Rouwenhorst \n',znums(end),tF/tR)
fprintf('(%2.6f s against %2.6f s). Whether that is worth paying is what the accuracy figures answer, \n',tF,tR)
fprintf('and the answer differs by calibration: Farmer-Toda wins at moderate persistence and \n')
fprintf('Rouwenhorst takes over as rho approaches one. \n')

if any(outputP2.Tauchen.moderate.nrepsused==1)
    fprintf('\nNote: some configs exceeded the warm-up time threshold and were timed once rather than %i \n',calib.nreps)
    fprintf('      times, so those cells are noisier. \n')
end
fprintf('\nNo configs were skipped for size in this block (that only arises in the multi-dimensional \n')
fprintf('blocks, where prod(znum) drives the cost). \n')

end
