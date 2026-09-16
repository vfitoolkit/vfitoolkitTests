function output=DiscP0_MarkovChainMoments_vs_panel(calibP0,figure_c)
% P0: validate MarkovChainMoments() against simulated panel data
%
% This is the block that earns the bank the right to use MarkovChainMoments as its measuring
% instrument everywhere else.
%
% The trick is n_a=1. With a single asset point the household has nothing to choose, so the model
% degenerates to "an agent carried along by the exogenous shock" and a simulated panel of z is a
% PURE MONTE CARLO DRAW from the discretized chain - produced by an entirely separate code path
% (SimPanelValues walks pi_z with rand) from the one MarkovChainMoments uses (an eigenvector solve).
% Two independent routes to the same population quantity.
%
% These are Monte Carlo checks and are never exact:
%   - SimPanelValues is not reproducible run to run (the parfor workers do not take the client's
%     rng), so two runs of this subcode will not agree to the last digit. Nothing here may be
%     written as an exact zero.
%   - The tolerance is COMPUTED, not chosen: each check is |panel - chain| < kse*se, with se
%     estimated from the panel itself. Every check prints its achieved se, so the diary records
%     what the block was actually able to detect.
%   - kse=4 is about a 1-in-16,000 false positive rate per check. With ~40 checks here that is
%     roughly a 0.25% chance of one spurious red per run, so a LONE red here means re-run once
%     before believing it.

fprintf('\n========== P0: MarkovChainMoments against simulated panel data ========== \n')
fprintf('Monte Carlo checks: tolerance is kse=%i standard errors, se estimated from the panel. \n',calibP0.kse)
fprintf('numbersims=%i, simperiods=%i \n',calibP0.numbersims,calibP0.simperiods)

output=struct();
N=calibP0.numbersims;
kse=calibP0.kse;

%% The three representative chains
% (a) Farmer-Toda AR(1) at moderate persistence
[zg{1},pz{1}]=discretizeAR1_FarmerToda(calibP0.panel.mew,calibP0.panel.rho,calibP0.panel.sigma,calibP0.panel.znum);
nm{1}='FarmerToda AR(1), rho=0.7';
% (b) Rouwenhorst at high persistence (a different construction, and a slower-mixing chain)
[zg{2},pz{2}]=discretizeAR1_Rouwenhorst(calibP0.panel.mew,calibP0.panel.rho_persistent,calibP0.panel.sigma,calibP0.panel.znum);
nm{2}='Rouwenhorst AR(1), rho=0.95';
% (c) The hand-built two-state chain: deliberately asymmetric, so a symmetric bug cannot hide
p=calibP0.two.p; q=calibP0.two.q;
zg{3}=[calibP0.two.z1; calibP0.two.z2];
pz{3}=[1-p, p; q, 1-q];
nm{3}='two-state, asymmetric';

%% The trivial n_a=1 model
% Note: if the toolkit cannot solve n_a=1, that is itself a finding rather than a setup problem.
n_d=0; d_grid=[];
n_a=1; a_grid=1;
Params=calibP0.Params;
DiscountFactorParamNames={'beta'};
ReturnFn=@(aprime,a,z) z; % payoff is irrelevant here; there is nothing to choose
FnsToEvaluate.zvalue=@(aprime,a,z) z;

storemean=zeros(3,2); storesd=zeros(3,2); storeac=zeros(3,2);

%% Warm up the parallel pool BEFORE any timing
% SimPanelValues runs its draws on parfor workers, so the first call of the session pays for
% starting the pool. Without this, the first chain's panel timing is dominated by parpool startup
% and reads as a 45x difference from the others when the panels are in fact the same size.
fprintf('warming up the parallel pool (its startup would otherwise land in the first timed call) \n')
vfoptionsW=struct(); simoptionsW=struct();
z_gridW=zg{1}; pi_zW=pz{1}; n_zW=length(z_gridW);
[~,PolicyW]=ValueFnIter_InfHorz(n_d,n_a,n_zW,d_grid,a_grid,z_gridW,pi_zW,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsW);
StationaryDistW=StationaryDist_InfHorz(PolicyW,n_d,n_a,n_zW,pi_zW,simoptionsW,Params,[]);
simoptionsW.numbersims=100; simoptionsW.simperiods=5; simoptionsW.burnin=0;
SimPanelValues_InfHorz(StationaryDistW,PolicyW,FnsToEvaluate,[],Params,n_d,n_a,n_zW,d_grid,a_grid,z_gridW,pi_zW,simoptionsW);
clear vfoptionsW simoptionsW z_gridW pi_zW n_zW PolicyW StationaryDistW

for c_c=1:3
    z_grid=zg{c_c}; pi_z=pz{c_c}; n_z=length(z_grid);
    fprintf('\n--- %s (znum=%i) --- \n',nm{c_c},n_z)

    % What the instrument says
    tic;
    [mcmean,mcvar,mcautocorr,mcstatdist]=MarkovChainMoments(z_grid,pi_z);
    t_instrument=toc;

    % What a simulation says
    vfoptions=struct(); simoptions=struct();
    [~,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params,[]);
    simoptionsPanel=simoptions;
    simoptionsPanel.numbersims=N;
    simoptionsPanel.simperiods=calibP0.simperiods;
    simoptionsPanel.burnin=0; % InitialDist is the stationary dist, so there is nothing to burn in
    tic;
    SimPanel=SimPanelValues_InfHorz(StationaryDist,Policy,FnsToEvaluate,[],Params,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,simoptionsPanel);
    t_panel=toc;

    % The period-1 cross-section is a straight draw from the stationary distribution, so its
    % observations are INDEPENDENT across agents and N_effective is exactly numbersims. Using the
    % pooled panel instead would inflate the apparent sample size, because observations along one
    % path are serially correlated.
    x1=SimPanel.zvalue(1,:);
    x1=x1(:);

    % mean
    se_mean=std(x1)/sqrt(N);
    d_mean=abs(mean(x1)-mcmean);
    fprintf('mean: chain %2.6f, panel %2.6f, |diff| %2.3e, se %2.3e [MC, pass if |diff|<%i*se: %i] \n',mcmean,mean(x1),d_mean,se_mean,kse,d_mean<kse*se_mean)

    % standard deviation. se(sigmahat)=sigma/sqrt(2N)
    se_sd=std(x1)/sqrt(2*N);
    d_sd=abs(std(x1,1)-sqrt(mcvar));
    fprintf('std dev: chain %2.6f, panel %2.6f, |diff| %2.3e, se %2.3e [MC, pass if |diff|<%i*se: %i] \n',sqrt(mcvar),std(x1,1),d_sd,se_sd,kse,d_sd<kse*se_sd)

    % autocorrelation, from the within-agent lag-1 pairs pooled over agents.
    % N_effective is taken as numbersims (not numbersims*(simperiods-1)): pairs within one path
    % are not independent, so this is the conservative count.
    xlag=SimPanel.zvalue(1:end-1,:); xnext=SimPanel.zvalue(2:end,:);
    temp=corrcoef(xlag(:),xnext(:)); r=temp(1,2);
    se_ac=(1-r^2)/sqrt(N);
    d_ac=abs(r-mcautocorr);
    fprintf('autocorrelation: chain %2.6f, panel %2.6f, |diff| %2.3e, se %2.3e [MC, pass if |diff|<%i*se: %i] \n',mcautocorr,r,d_ac,se_ac,kse,d_ac<kse*se_ac)
    % This is the check that catches a transposed pi_z. The stationary distribution is unchanged
    % under transposition whenever pi_z is (near-)reversible, so mean and std dev can both pass
    % while the dynamics run backwards.

    % the stationary distribution itself, point by point. se for a proportion is sqrt(p(1-p)/N).
    % This is what localises WHICH grid point is wrong when one of the above fails.
    empfreq=zeros(n_z,1);
    for z_c=1:n_z
        empfreq(z_c)=sum(abs(x1-z_grid(z_c))<10^(-10))/N;
    end
    se_freq=sqrt(mcstatdist(:).*(1-mcstatdist(:))/N);
    d_freq=abs(empfreq-mcstatdist(:));
    fprintf('stationary dist, worst grid point: |diff| %2.3e vs %i*se %2.3e [MC, pass: %i] \n',max(d_freq),kse,max(kse*se_freq),all(d_freq<kse*se_freq))
    fprintf('   (empirical frequencies sum to %2.6f; if this is not 1, the panel left the grid) \n',sum(empfreq))

    storemean(c_c,:)=[mcmean,mean(x1)];
    storesd(c_c,:)=[sqrt(mcvar),std(x1,1)];
    storeac(c_c,:)=[mcautocorr,r];
    output.chain(c_c).name=nm{c_c};
    output.chain(c_c).mean=[mcmean,mean(x1),se_mean];
    output.chain(c_c).sd=[sqrt(mcvar),std(x1,1),se_sd];
    output.chain(c_c).autocorrelation=[mcautocorr,r,se_ac];
    output.chain(c_c).statdist=[mcstatdist(:),empfreq(:),se_freq(:)];
    output.chain(c_c).runtime=[t_instrument,t_panel]; % seconds: MarkovChainMoments, SimPanelValues
end

%% Figure
figure(figure_c)
subplot(1,3,1); bar(storemean); title('mean'); set(gca,'XTickLabel',{'FT','Rouw','2-state'})
legend('MarkovChainMoments','panel','Location','best')
subplot(1,3,2); bar(storesd); title('std dev'); set(gca,'XTickLabel',{'FT','Rouw','2-state'})
subplot(1,3,3); bar(storeac); title('autocorrelation'); set(gca,'XTickLabel',{'FT','Rouw','2-state'})

end
