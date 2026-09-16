function DiscP1_runtimetable(outputP1)
% Runtime table for P1 (iid normal)
%
% Runtimes are REPORTED, never asserted: they are machine-, MATLAB-version- and thermal-dependent,
% and a bank that goes red because a laptop throttled is a bank people stop running. What they are
% for is (i) telling a user which method to pick when accuracy is a tie - the spread across a block
% is often two orders of magnitude, which is decision-relevant - and (ii) making an accidental
% complexity regression visible to a human reading the diary.
%
% These are not extra runs: the timed calls are the same calls the accuracy sweep read (§5 of the
% proposal), so the runtime table and the accuracy figure share an x-axis by construction.

fprintf('\n========== Runtimes for discretizing an iid normal ========== \n')
v=version;
try
    g=gpuDevice; gname=g.Name;
catch
    gname='no gpu';
end
fprintf('(median of up to %i timed calls after 1 discarded warm-up; MATLAB %s; GPU: %s) \n',outputP1.IIDNormal_Tauchen.nrepsused(1),v,gname)

znums=outputP1.IIDNormal_Tauchen.znums;
fprintf('%-38s',' ');
for c_c=1:length(znums)
    fprintf('%11s',['znum=',num2str(znums(c_c))]);
end
fprintf(' \n');

% The rows. Written out with an explicit loop rather than a local subfunction, per the toolkit's
% no-helper-functions rule.
rowname={'discretizeIIDNormal_Tauchen','discretizeIIDNormal_TanakaToda','discretizeIID_Tauchen','discretizeIID_TanakaToda'};
rowtime=[outputP1.IIDNormal_Tauchen.runtime; outputP1.IIDNormal_TanakaToda.runtime; outputP1.IID_Tauchen.runtime; outputP1.IID_TanakaToda.runtime];
for m_c=1:4
    rowname{4+m_c}=['  ',outputP1.AR1methods_rho0.nmlist{m_c},' (rho=0)'];
    rowtime(4+m_c,:)=outputP1.AR1methods_rho0.runtime(m_c,:);
end
for m_c=1:length(rowname)
    fprintf('%-38s',rowname{m_c});
    for c_c=1:length(znums)
        fprintf('%11.6f',rowtime(m_c,c_c));
    end
    fprintf(' \n');
end

% Any config whose warm-up call exceeded the threshold was timed once rather than nreps times;
% flag those, because a single timing is noisier and the reader should know which cells are which.
if any(outputP1.IIDNormal_Tauchen.nrepsused==1) || any(outputP1.IIDNormal_TanakaToda.nrepsused==1)
    fprintf('\nNote: some configs exceeded the warm-up time threshold and were timed once rather than \n')
    fprintf('      %i times, so those cells are noisier. \n',outputP1.IIDNormal_Tauchen.nrepsused(1))
end
fprintf('\nNo configs were skipped for size in this block (that only arises in the multi-dimensional \n')
fprintf('blocks, where prod(znum) drives the cost). \n')

end
