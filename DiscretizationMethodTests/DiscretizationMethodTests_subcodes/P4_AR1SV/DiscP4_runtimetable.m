function DiscP4_runtimetable(outputP4,calib,znums)
% Runtime table for P4 (AR(1) with stochastic volatility)
%
% Reported, never asserted. This is the first block whose table is TWO-DIMENSIONAL, because the
% sweep is nested over xnum and znum, and the first where configs are skipped for size - so the
% skipped cells are marked rather than left blank. A nested sweep that quietly drops cells looks
% identical in a diary to one that ran them all.

fprintf('\n========== Runtimes for discretizing an AR(1) with stochastic volatility ========== \n')
v=version;
try
    g=gpuDevice; gname=g.Name;
catch
    gname='no gpu';
end
fprintf('(median of up to %i timed calls after 1 discarded warm-up; MATLAB %s; GPU: %s) \n',calib.nreps,v,gname)
fprintf('cells with xnum*znum above the ceiling of %i are skipped and shown as "skip" \n',calib.Ncap)

nmlist={'discretizeAR1wSV_FarmerToda','discretizeAR1wSV_Tauchen'};
fldlist={'AR1wSV_FarmerToda','AR1wSV_Tauchen'};
nz=length(znums);
for m_c=1:2
    fprintf('\n--- %s (rows xnum, columns znum) --- \n',nmlist{m_c})
    fprintf('%8s',' ');
    for c_c=1:nz
        fprintf('%11i',znums(c_c));
    end
    fprintf(' \n');
    t=outputP4.(fldlist{m_c}).runtime;
    sk=outputP4.(fldlist{m_c}).skipped;
    for x_c=1:nz
        fprintf('%8i',znums(x_c));
        for z_c=1:nz
            if sk(x_c,z_c)
                fprintf('%11s','skip');
            else
                fprintf('%11.6f',t(x_c,z_c));
            end
        end
        fprintf(' \n');
    end
end

nskip=sum(outputP4.AR1wSV_FarmerToda.skipped(:));
fprintf('\n%i of %i cells were skipped for size. The cost here is driven by the joint transition \n',nskip,nz*nz)
fprintf('matrix, which is (xnum*znum)-by-(xnum*znum), so it grows as the SQUARE of the product of \n')
fprintf('the two grid sizes rather than of either one. \n')

end
