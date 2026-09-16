function DiscP3_runtimetable(outputP3,calib,znums)
% Runtime table for P3 (AR(1) with gaussian-mixture innovations)
%
% Reported, never asserted. Measured on the same calls the accuracy sweep read.

fprintf('\n========== Runtimes for discretizing an AR(1) with gaussian-mixture innovations ========== \n')
v=version;
try
    g=gpuDevice; gname=g.Name;
catch
    gname='no gpu';
end
fprintf('(median of up to %i timed calls after 1 discarded warm-up; MATLAB %s; GPU: %s) \n',calib.nreps,v,gname)

fprintf('%-34s',' ');
for c_c=1:length(znums)
    fprintf('%11s',['znum=',num2str(znums(c_c))]);
end
fprintf(' \n');

rowname={'discretizeAR1wGM_FarmerToda','discretizeAR1wGM_Tauchen'};
rowtime=[outputP3.AR1wGM_FarmerToda.runtime; outputP3.AR1wGM_Tauchen.runtime];
for m_c=1:2
    fprintf('%-34s',rowname{m_c});
    for c_c=1:length(znums)
        fprintf('%11.6f',rowtime(m_c,c_c));
    end
    fprintf(' \n');
end
% and the gaussian comparators, which are what a user would otherwise reach for
gnm=outputP3.gaussianmethods.nmlist;
for m_c=1:4
    fprintf('%-34s',['  ',gnm{m_c},' (misspecified)']);
    for c_c=1:length(znums)
        fprintf('%11.6f',outputP3.gaussianmethods.runtime(m_c,c_c));
    end
    fprintf(' \n');
end

fprintf('\nThe comparison to read off this table together with the accuracy figures: the mixture \n')
fprintf('methods cost more, and what that buys is the skewness and excess kurtosis, which no amount \n')
fprintf('of grid refinement will get from a gaussian method. If a model does not care about the \n')
fprintf('higher moments, the table says the mixture machinery is not worth paying for. \n')

end
