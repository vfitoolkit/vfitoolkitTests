function output=DiscP0_MarkovChainMoments_FHorz(calibP0,figure_c)
% P0: closed-form validation of MarkovChainMoments_FHorz()
%
% Built on a HAND-MADE age-dependent chain rather than on the output of a discretization command,
% for two reasons: the chain's age profiles are then known exactly rather than approximately, and
% a failure is unambiguously the instrument.
%
% What this does and does not test. The moments are recomputed here with explicit loops, using the
% same definitions the instrument uses, so this checks the instrument's BOOKKEEPING - which age's
% grid goes with which slice of pi_z_J, whether the age indexing is off by one, whether the lag
% pairing in the autocorrelation is the right way round - rather than the definitions themselves.
% That is the failure mode worth guarding: with a shifted age index every profile still looks
% plausible.
%
% The grid VARIES with age here on purpose. With one grid for all ages, using the wrong age's grid
% is invisible.

fprintf('\n========== P0: MarkovChainMoments_FHorz against a hand-built chain ========== \n')

output=struct();

znum=calibP0.agedep.znum;
N_j=calibP0.agedep.N_j;
jequaloneDistz=calibP0.agedep.jequaloneDistz;

%% Build the chain by hand
% Grid: the base grid, shifted and stretched by age, so every age is distinguishable
z_grid_J=zeros(znum,N_j);
for jj=1:N_j
    z_grid_J(:,jj)=calibP0.agedep.z_grid*(1+0.1*jj)+0.05*jj;
end

% Transition matrices: pi_z_J(:,:,jj) is from age jj to age jj+1, so there are N_j-1 of them.
% Each age gets a different matrix, again so that a wrong slice is visible.
pi_z_J=zeros(znum,znum,N_j-1);
for jj=1:N_j-1
    w=0.5+0.05*jj; % how much mass stays put at age jj, varying with age
    P=zeros(znum,znum);
    for z_c=1:znum
        P(z_c,:)=(1-w)/(znum-1);
        P(z_c,z_c)=w;
    end
    % Deliberately make it non-symmetric, so that a transposed slice is detectable
    P(1,znum)=P(1,znum)+0.1; P(1,1)=P(1,1)-0.1;
    pi_z_J(:,:,jj)=P;
end
% Sanity: rows must be distributions, otherwise this subcode is testing nothing
fprintf('hand-built chain: rows of pi_z_J sum to one [T1], this should be zero: %2.8e \n',max(abs(sum(pi_z_J,2)-1),[],'all'))
fprintf('hand-built chain: jequaloneDistz sums to one [T1], this should be zero: %2.8e \n',abs(sum(jequaloneDistz)-1))

%% Recompute the age profiles here, with explicit loops
statdist_true=zeros(znum,N_j);
mean_true=zeros(1,N_j);
var_true=zeros(1,N_j);
autocorr_true=nan(1,N_j); % j=1 has no lag, so it stays nan

statdist_true(:,1)=jequaloneDistz;
mean_true(1)=sum(statdist_true(:,1).*z_grid_J(:,1));
var_true(1)=sum(statdist_true(:,1).*(z_grid_J(:,1)-mean_true(1)).^2);
for jj=2:N_j
    % Push the distribution forward one age
    for z_c=1:znum
        statdist_true(z_c,jj)=sum(statdist_true(:,jj-1).*pi_z_J(:,z_c,jj-1));
    end
    mean_true(jj)=sum(statdist_true(:,jj).*z_grid_J(:,jj));
    var_true(jj)=sum(statdist_true(:,jj).*(z_grid_J(:,jj)-mean_true(jj)).^2);
    % Covariance between age jj-1 and age jj, summed over the joint distribution of the pair
    covar=0;
    for z_c=1:znum       % state at age jj-1
        for zprime_c=1:znum  % state at age jj
            jointprob=statdist_true(z_c,jj-1)*pi_z_J(z_c,zprime_c,jj-1);
            covar=covar+jointprob*(z_grid_J(z_c,jj-1)-mean_true(jj-1))*(z_grid_J(zprime_c,jj)-mean_true(jj));
        end
    end
    autocorr_true(jj)=covar/(sqrt(var_true(jj-1))*sqrt(var_true(jj)));
end

%% Compare
[mcmean,mcvar,mcautocorr,mcstatdist]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);

fprintf('age profile of the stationary dist [T1], this should be zero: %2.8e \n',max(abs(mcstatdist(:)-statdist_true(:))))
fprintf('age profile of the mean [T1], this should be zero: %2.8e \n',max(abs(mcmean(:)-mean_true(:))))
fprintf('age profile of the variance [T1], this should be zero: %2.8e \n',max(abs(mcvar(:)-var_true(:))))
fprintf('age profile of the autocorrelation (ages 2:N_j) [T1], this should be zero: %2.8e \n',max(abs(mcautocorr(2:end)-autocorr_true(2:end))))
output.statdist=[mcstatdist(:),statdist_true(:)];
output.mean=[mcmean(:),mean_true(:)];
output.variance=[mcvar(:),var_true(:)];
output.autocorrelation=[mcautocorr(:),autocorr_true(:)];

% Age 1 must reproduce the initial distribution exactly; this is documented behaviour
fprintf('age 1 statdist reproduces jequaloneDistz [T0], this should be zero: %2.8e \n',max(abs(mcstatdist(:,1)-jequaloneDistz(:))))

%% The age index must actually be used: shifting the chain by one age must change the answer
% If this prints zero, the checks above are not testing which slice goes with which age.
pi_z_J_shifted=pi_z_J(:,:,[2:end,end]);
[~,~,~,mcstatdist_shift]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J_shifted,jequaloneDistz);
fprintf('shifting pi_z_J by one age must CHANGE the profile, this should NOT be zero: %2.8e \n',max(abs(mcstatdist_shift(:)-mcstatdist(:))))

%% pi_z_J with N_j-1 slices (the current convention) and with N_j slices (the old one, still
%% accepted) must give the same answer: only slices 1..N_j-1 are ever read.
pi_z_J_oldshape=zeros(znum,znum,N_j);
pi_z_J_oldshape(:,:,1:N_j-1)=pi_z_J;
pi_z_J_oldshape(:,:,N_j)=ones(znum,znum)/znum; % the meaningless padding slice of the old convention
[m2,v2,a2,s2]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J_oldshape,jequaloneDistz);
fprintf('pi_z_J with N_j-1 vs N_j slices: statdist [T0], this should be zero: %2.8e \n',max(abs(s2(:)-mcstatdist(:))))
fprintf('pi_z_J with N_j-1 vs N_j slices: mean [T0], this should be zero: %2.8e \n',max(abs(m2(:)-mcmean(:))))
fprintf('pi_z_J with N_j-1 vs N_j slices: variance [T0], this should be zero: %2.8e \n',max(abs(v2(:)-mcvar(:))))
fprintf('pi_z_J with N_j-1 vs N_j slices: autocorrelation [T0], this should be zero: %2.8e \n',max(abs(a2(2:end)-mcautocorr(2:end))))

%% Option sweep
% simoptions.z_grid_J / pi_z_J / jequaloneDist overwrite the positional inputs; passing the same
% objects both ways must be a no-op
simoptions=struct();
simoptions.z_grid_J=z_grid_J;
simoptions.pi_z_J=pi_z_J;
simoptions.jequaloneDist=jequaloneDistz;
[m3,v3,a3,s3]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz,simoptions);
fprintf('simoptions overwrite with the same objects: statdist [T0], this should be zero: %2.8e \n',max(abs(s3(:)-mcstatdist(:))))
fprintf('simoptions overwrite with the same objects: mean [T0], this should be zero: %2.8e \n',max(abs(m3(:)-mcmean(:))))
fprintf('simoptions overwrite with the same objects: variance [T0], this should be zero: %2.8e \n',max(abs(v3(:)-mcvar(:))))
fprintf('simoptions overwrite with the same objects: autocorrelation [T0], this should be zero: %2.8e \n',max(abs(a3(2:end)-mcautocorr(2:end))))

% calcautocorrelation=0 must not disturb the other three
mcopts=struct(); mcopts.calcautocorrelation=0;
[m4,v4,~,s4]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz,struct(),mcopts);
fprintf('calcautocorrelation=0: mean unchanged [T0], this should be zero: %2.8e \n',max(abs(m4(:)-mcmean(:))))
fprintf('calcautocorrelation=0: variance unchanged [T0], this should be zero: %2.8e \n',max(abs(v4(:)-mcvar(:))))
fprintf('calcautocorrelation=0: statdist unchanged [T0], this should be zero: %2.8e \n',max(abs(s4(:)-mcstatdist(:))))

%% Figure
figure(figure_c)
subplot(2,2,1); plot(1:N_j,mean_true,'k-',1:N_j,mcmean,'ro')
title('mean by age'); xlabel('age j'); legend('recomputed here','MarkovChainMoments\_FHorz','Location','best')
subplot(2,2,2); plot(1:N_j,var_true,'k-',1:N_j,mcvar,'ro')
title('variance by age'); xlabel('age j')
subplot(2,2,3); plot(2:N_j,autocorr_true(2:end),'k-',2:N_j,mcautocorr(2:end),'ro')
title('autocorrelation (j-1,j)'); xlabel('age j')
subplot(2,2,4); plot(1:N_j,abs(mcmean-mean_true),'o-',1:N_j,abs(mcvar-var_true),'s-')
set(gca,'YScale','log'); title('|error| by age'); xlabel('age j'); legend('mean','variance','Location','best')

end
