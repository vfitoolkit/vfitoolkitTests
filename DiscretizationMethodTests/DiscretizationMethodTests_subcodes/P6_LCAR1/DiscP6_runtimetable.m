function DiscP6_runtimetable(outputP6,calib,znums)
% Runtime table for P6 (life-cycle AR(1) with normal innovations)
%
% Reported, never asserted. J is fixed at the calibration's value, so these numbers scale roughly
% linearly in J: every one of these commands does per-age work and nothing that couples the ages.
% What differs between them is the per-age cost - KFTT solves a maximum entropy problem for every
% grid point at every age, the two Fella-Gallipoli-Pan methods do closed-form arithmetic.

fprintf('\n========== Runtimes for discretizing a life-cycle AR(1) with normal innovations ========== \n')
v=version;
try
    g=gpuDevice; gname=g.Name;
catch
    gname='no gpu';
end
fprintf('(median of up to %i timed calls after 1 discarded warm-up; J=%i; MATLAB %s; GPU: %s) \n',calib.nreps,calib.J,v,gname)

nmlist={'discretizeLifeCycleAR1_KFTT','discretizeLifeCycleAR1_FellaGallipoliPan','discretizeLifeCycleAR1_FellaGallipoliPanTauchen'};
fldlist={'LCAR1_KFTT','LCAR1_FGP','LCAR1_FGPTauchen'};
nz=length(znums);
fprintf('\n%10s','znum');
for m_c=1:3
    fprintf('%16s',fldlist{m_c});
end
fprintf('%14s \n','KFTT/FGP');
for c_c=1:nz
    fprintf('%10i',znums(c_c));
    t=zeros(1,3);
    for m_c=1:3
        t(m_c)=outputP6.(fldlist{m_c}).runtime(c_c);
        fprintf('%16.6f',t(m_c));
    end
    fprintf('%14.1f \n',t(1)/t(2));
end

fprintf('\nfull command names: \n');
for m_c=1:3
    fprintf('   %-18s %s \n',fldlist{m_c},nmlist{m_c});
end
fprintf('\nThe last column is the price of the KFTT method: it solves a maximum entropy problem for \n')
fprintf('every grid point at every age, where the two Fella-Gallipoli-Pan methods evaluate a closed \n')
fprintf('form. Whether that price is worth paying is the comparison in DiscP6_compare, and the two \n')
fprintf('FGP commands themselves print, on every call, that it is. \n')

end
