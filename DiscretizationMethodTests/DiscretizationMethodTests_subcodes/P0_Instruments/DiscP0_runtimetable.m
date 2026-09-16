function DiscP0_runtimetable(outputP0)
% Runtime table for P0.
%
% Runtimes in this bank are REPORTED, never asserted: they are machine-, MATLAB-version- and
% thermal-dependent, and a bank that goes red because a laptop throttled is a bank people stop
% running. What they are for is (i) telling a user which method to pick when accuracy is a tie,
% and (ii) making an accidental complexity regression visible to a human reading the diary.
%
% P0's table is short but worth having: the panel validation is by far the slowest thing in the
% bank and it runs first, so this is where "why is the bank slow to start" gets answered.

fprintf('\n========== Runtimes for P0 (validating the instruments) ========== \n')
v=version;
try
    g=gpuDevice; gname=g.Name;
catch
    gname='no gpu';
end
fprintf('(MATLAB %s; GPU: %s) \n',v,gname)
fprintf('%-46s %12s %12s \n','','instrument','panel')

for c_c=1:length(outputP0.vs_panel.chain)
    fprintf('%-46s %12.5f %12.5f \n',outputP0.vs_panel.chain(c_c).name,outputP0.vs_panel.chain(c_c).runtime(1),outputP0.vs_panel.chain(c_c).runtime(2))
end
fprintf('%-46s %12.5f %12.5f \n','life-cycle KFTT (FHorz)',outputP0.vs_panel_FHorz.runtime(1),outputP0.vs_panel_FHorz.runtime(2))

fprintf('\nNote: the "instrument" column is MarkovChainMoments(_FHorz), the "panel" column is \n')
fprintf('      SimPanelValues(_FHorz_Case1). The panel column is what makes the bank slow to start; \n')
fprintf('      it is the price of validating the instrument before anything else depends on it. \n')

end
