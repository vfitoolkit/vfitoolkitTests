function DiscP7_runtimetable(outputP7,calib,znums)
% Runtime table for P7 (life-cycle AR(1) with gaussian-mixture innovations)
%
% Reported, never asserted. This is the most expensive block in the bank, and the reason is
% structural rather than incidental: discretizeLifeCycleAR1wGM_KFTT solves a maximum entropy problem
% for every grid point at every age, so the work is znum*J solves per call - 2091 of them at
% znum=51, J=41 - where discretizeLifeCycleAR1wGM_Tauchen evaluates a closed-form mixture cdf and
% carries no solver at all. The ratio in the last column is what the KFTT method costs, and whether
% it is worth paying is the comparison in DiscP7_compare.

fprintf('\n========== Runtimes for discretizing a life-cycle AR(1) with gaussian-mixture innovations ========== \n')
v=version;
try
    g=gpuDevice; gname=g.Name;
catch
    gname='no gpu';
end
fprintf('(median of up to %i timed calls after 1 discarded warm-up; J=%i; MATLAB %s; GPU: %s) \n',calib.nreps,calib.J,v,gname)

nz=length(znums);
fprintf('\n%10s%18s%18s%14s%16s \n','znum','wGM_KFTT','wGM_Tauchen','KFTT/Tauchen','entropy solves');
for c_c=1:nz
    tK=outputP7.LCAR1wGM_KFTT.runtime(c_c);
    tT=outputP7.LCAR1wGM_Tauchen.runtime(c_c);
    fprintf('%10i%18.6f%18.6f%14.1f%16i \n',znums(c_c),tK,tT,tK/tT,znums(c_c)*calib.J);
end

fprintf('\nThe last column is the number of maximum entropy problems the KFTT method solves per call. \n')
fprintf('It is linear in both znum and J, so this block scales with the product - which is why P7 \n')
fprintf('is the slowest in the bank and why its sweeps stop at znum=51 rather than 101. \n')

end
