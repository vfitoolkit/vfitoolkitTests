function output=DiscP0_MarkovChainMoments_FHorz_vs_panel(calibP0,figure_c)
% P0: validate MarkovChainMoments_FHorz() against simulated panel data
%
% Same idea as the InfHorz half: n_a=1, so a simulated panel of z is a pure Monte Carlo draw from
% the discretized age-dependent chain, produced by a code path entirely separate from the one the
% instrument uses.
%
% What this half adds over the InfHorz one is the AGE INDEXING. If pi_z_J(:,:,j) were being read
% as the transition INTO age j rather than OUT OF it, every age profile would be shifted by one
% and nothing else in this bank would notice. It also exercises jequaloneDistz, which is otherwise
% only ever checked for summing to one.
%
% Monte Carlo throughout; see the InfHorz subcode's header for the tolerance rules. Because there
% are N_j ages, the per-age lines would swamp the diary, so what is printed is the WORST age for
% each moment (in units of se) and the age at which it occurred. The full age profiles are in the
% figure and in the returned struct.

fprintf('\n========== P0: MarkovChainMoments_FHorz against simulated panel data ========== \n')
fprintf('Monte Carlo checks: tolerance is kse=%i standard errors, se estimated from the panel. \n',calibP0.kse)

output=struct();
N=calibP0.numbersims;
kse=calibP0.kse;
N_j=calibP0.panelJ.N_j;
znum=calibP0.panelJ.znum;

%% The life-cycle chain
kfttoptions=struct();
[z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleAR1_KFTT(calibP0.panelJ.mew,calibP0.panelJ.rho,calibP0.panelJ.sigma,znum,N_j,kfttoptions);
fprintf('chain built: size(pi_z_J,3)=%i, expected N_j-1=%i \n',size(pi_z_J,3),N_j-1)

%% What the instrument says
tic;
[mcmean,mcvar,mcautocorr,mcstatdist]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);
t_instrument=toc;

%% What a simulation says
n_d=0; d_grid=[];
n_a=1; a_grid=1;
n_z=znum;
Params=calibP0.Params;
DiscountFactorParamNames={'beta'};
ReturnFn=@(aprime,a,z) z; % payoff is irrelevant; there is nothing to choose
FnsToEvaluate.zvalue=@(aprime,a,z) z;

vfoptions=struct(); simoptions=struct();
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);

jequaloneDist=zeros(n_a,n_z);
jequaloneDist(1,:)=jequaloneDistz(:)'; % no assets to speak of, so the initial dist is just the z one

simoptionsPanel=simoptions;
simoptionsPanel.numbersims=N;
simoptionsPanel.simperiods=N_j;
tic;
SimPanel=SimPanelValues_FHorz_Case1(jequaloneDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,pi_z_J,simoptionsPanel);
t_panel=toc;

zpanel=SimPanel.zvalue; % [N_j, numbersims]

%% Age-conditional mean and standard deviation
% Agents are independent, so at each age the cross-section is N independent draws and N_effective
% is exactly numbersims.
panelmean=zeros(1,N_j); panelsd=zeros(1,N_j);
se_mean=zeros(1,N_j); se_sd=zeros(1,N_j);
for jj=1:N_j
    xj=zpanel(jj,:); xj=xj(:);
    panelmean(jj)=mean(xj);
    panelsd(jj)=std(xj,1);
    se_mean(jj)=std(xj)/sqrt(N);
    se_sd(jj)=std(xj)/sqrt(2*N);
end
d_mean=abs(panelmean-mcmean(:)');
d_sd=abs(panelsd-sqrt(mcvar(:)'));
ratio_mean=d_mean./se_mean;
ratio_sd=d_sd./se_sd;
[worstm,jm]=max(ratio_mean);
[worsts,js]=max(ratio_sd);
fprintf('age-conditional mean: worst age is j=%i at %2.2f se [MC, pass if <%i: %i] \n',jm,worstm,kse,worstm<kse)
fprintf('   at that age: chain %2.6f, panel %2.6f, se %2.3e \n',mcmean(jm),panelmean(jm),se_mean(jm))
fprintf('age-conditional std dev: worst age is j=%i at %2.2f se [MC, pass if <%i: %i] \n',js,worsts,kse,worsts<kse)
fprintf('   at that age: chain %2.6f, panel %2.6f, se %2.3e \n',sqrt(mcvar(js)),panelsd(js),se_sd(js))

%% Age-to-age autocorrelation. This is the one that catches a shifted age index.
panelac=nan(1,N_j); se_ac=nan(1,N_j);
for jj=2:N_j
    temp=corrcoef(zpanel(jj-1,:)',zpanel(jj,:)');
    panelac(jj)=temp(1,2);
    se_ac(jj)=(1-panelac(jj)^2)/sqrt(N);
end
d_ac=abs(panelac-mcautocorr(:)');
ratio_ac=d_ac./se_ac;
[worsta,ja]=max(ratio_ac(2:end)); ja=ja+1;
fprintf('age-to-age autocorrelation: worst age is j=%i at %2.2f se [MC, pass if <%i: %i] \n',ja,worsta,kse,worsta<kse)
fprintf('   at that age: chain %2.6f, panel %2.6f, se %2.3e \n',mcautocorr(ja),panelac(ja),se_ac(ja))

%% The age-1 cross-section is a straight draw from jequaloneDistz, point by point
empfreq1=zeros(znum,1);
x1=zpanel(1,:);
for z_c=1:znum
    empfreq1(z_c)=sum(abs(x1-z_grid_J(z_c,1))<10^(-10))/N;
end
se_freq=sqrt(jequaloneDistz(:).*(1-jequaloneDistz(:))/N);
d_freq=abs(empfreq1-jequaloneDistz(:));
fprintf('age 1 cross-section vs jequaloneDistz, worst grid point: |diff| %2.3e vs %i*se %2.3e [MC, pass: %i] \n',max(d_freq),kse,max(kse*se_freq),all(d_freq<kse*se_freq))
fprintf('   (empirical frequencies sum to %2.6f; if this is not 1, the panel left the grid) \n',sum(empfreq1))

output.mean=[mcmean(:),panelmean(:),se_mean(:)];
output.sd=[sqrt(mcvar(:)),panelsd(:),se_sd(:)];
output.autocorrelation=[mcautocorr(:),panelac(:),se_ac(:)];
output.jequalonecheck=[jequaloneDistz(:),empfreq1(:),se_freq(:)];
output.statdist=mcstatdist;
output.runtime=[t_instrument,t_panel]; % seconds: MarkovChainMoments_FHorz, SimPanelValues_FHorz_Case1

%% Figure
figure(figure_c)
subplot(2,2,1); plot(1:N_j,mcmean,'k-',1:N_j,panelmean,'ro')
title('mean by age'); xlabel('age j'); legend('MarkovChainMoments\_FHorz','panel','Location','best')
subplot(2,2,2); plot(1:N_j,sqrt(mcvar),'k-',1:N_j,panelsd,'ro')
title('std dev by age'); xlabel('age j')
subplot(2,2,3); plot(2:N_j,mcautocorr(2:end),'k-',2:N_j,panelac(2:end),'ro')
title('autocorrelation (j-1,j)'); xlabel('age j')
subplot(2,2,4); plot(1:N_j,ratio_mean,'o-',1:N_j,ratio_sd,'s-',2:N_j,ratio_ac(2:end),'^-')
hold on; plot([1,N_j],[kse,kse],'k--'); hold off
title('|difference| in standard errors'); xlabel('age j'); legend('mean','std dev','autocorr','threshold','Location','best')

end
