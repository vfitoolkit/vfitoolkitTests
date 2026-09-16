function output=DiscP1_AR1methods_rho0(calib,znums,figure_c)
% P1: the four stationary AR(1) methods at rho=0
%
% An AR(1) with rho=0 IS an iid normal, so all four AR(1) methods discretize this block's process
% too. They are not swept here - P2 owns their option sweeps - they are run at defaults as
% comparison points, plus one check that is specific to rho=0 and sharp:
%
%   at rho=0 the conditional distribution does not depend on the lag, so EVERY ROW of pi_z must be
%   identical. That is exact for all four methods (Tauchen and Farmer-Toda because the conditional
%   mean is mew for every lag point, Rouwenhorst because p=q=1/2 gives the binomial matrix,
%   Tauchen-Hussey because EZprime is mew for every i). A method whose rows differ at rho=0 has
%   its lag indexing wrong, and nothing else in this bank would catch that.
%
% Note on mew: after the convention fix, all four take mew as the INTERCEPT, so E(z)=mew/(1-rho),
% which at rho=0 is just mew. The convention only bites at rho~=0, so pinning it is P2's job, not
% this subcode's.

fprintf('\n========== P1: the four AR(1) methods at rho=0 ========== \n')

output=struct();
mew=calib.mew; sigma=calib.sigma;
rho=0;
Tauchen_q=3;
nz=length(znums);

nmlist={'AR1_Tauchen','AR1_Rouwenhorst','AR1_FarmerToda','AR1_TauchenHussey'};
% Rouwenhorst matches the moments algebraically, so it is exact to machine precision. Farmer-Toda
% reaches them through an entropy solve, so "exactly" there means "to the solver's tolerance",
% which degrades with znum (measured: 1.3e-09 by znum=101). Hence two different thresholds.
exacttol=10^(-12);
entropytol=10^(-7);
err_mean=zeros(4,nz); err_var=zeros(4,nz); runtime=zeros(4,nz);

for c_c=1:nz
    znum=znums(c_c);
    for m_c=1:4
        % --- timing: one warm-up call, discarded, then nreps timed calls, report the median.
        % nreps drops to 1 for any config whose warm-up exceeded calib.timethreshold seconds.
        tic;
        if m_c==1
            [z_grid,pi_z]=discretizeAR1_Tauchen(mew,rho,sigma,znum,Tauchen_q,struct());
        elseif m_c==2
            [z_grid,pi_z]=discretizeAR1_Rouwenhorst(mew,rho,sigma,znum,struct());
        elseif m_c==3
            [z_grid,pi_z]=discretizeAR1_FarmerToda(mew,rho,sigma,znum,struct());
        else
            [z_grid,pi_z]=discretizeAR1_TauchenHussey(mew,rho,sigma,znum,struct());
        end
        twarm=toc;
        if twarm>calib.timethreshold
            nreps=1;
        else
            nreps=calib.nreps;
        end
        treps=zeros(1,nreps);
        for r_c=1:nreps
            tic;
            if m_c==1
                [z_grid,pi_z]=discretizeAR1_Tauchen(mew,rho,sigma,znum,Tauchen_q,struct());
            elseif m_c==2
                [z_grid,pi_z]=discretizeAR1_Rouwenhorst(mew,rho,sigma,znum,struct());
            elseif m_c==3
                [z_grid,pi_z]=discretizeAR1_FarmerToda(mew,rho,sigma,znum,struct());
            else
                [z_grid,pi_z]=discretizeAR1_TauchenHussey(mew,rho,sigma,znum,struct());
            end
            treps(r_c)=toc;
        end
        runtime(m_c,c_c)=median(treps);

        % --- invariants
        fprintf('%s, znum=%i: size of z_grid [T0], this should be zero: %i \n',nmlist{m_c},znum,any(size(z_grid)~=[znum,1]))
        fprintf('%s, znum=%i: size of pi_z [T0], this should be zero: %i \n',nmlist{m_c},znum,any(size(pi_z)~=[znum,znum]))
        fprintf('%s, znum=%i: z_grid is strictly ascending [T0], this should be one: %i \n',nmlist{m_c},znum,issorted(z_grid,'strictascend'))
        fprintf('%s, znum=%i: pi_z is in [0,1] [T0], this should be zero: %i \n',nmlist{m_c},znum,any(pi_z(:)<0)+any(pi_z(:)>1))
        fprintf('%s, znum=%i: rows of pi_z sum to one [T1], this should be zero: %2.8e \n',nmlist{m_c},znum,max(abs(sum(pi_z,2)-1)))
        fprintf('%s, znum=%i: no NaN or Inf [T0], this should be zero: %i \n',nmlist{m_c},znum,any(~isfinite(z_grid))+any(~isfinite(pi_z(:))))
        % symmetry: grid symmetric about the unconditional mean, pi_z centrosymmetric
        zstar=mew/(1-rho);
        fprintf('%s, znum=%i: z_grid symmetric about mew/(1-rho) [T1], this should be zero: %2.8e \n',nmlist{m_c},znum,max(abs((z_grid+flipud(z_grid))/2-zstar)))
        fprintf('%s, znum=%i: pi_z is centrosymmetric [T1], this should be zero: %2.8e \n',nmlist{m_c},znum,max(abs(pi_z-rot90(pi_z,2)),[],'all'))

        % --- the rho=0 check: every row identical
        fprintf('%s, znum=%i: at rho=0 every row of pi_z is identical [T1], this should be zero: %2.8e \n',nmlist{m_c},znum,max(abs(pi_z-repmat(pi_z(1,:),znum,1)),[],'all'))

        % --- accuracy. Rouwenhorst matches mean and variance exactly by construction, and
        % Farmer-Toda targets them, so those two are T1. Tauchen and Tauchen-Hussey approximate.
        pi_row=pi_z(1,:)';
        m=sum(pi_row.*z_grid);
        v=sum(pi_row.*(z_grid-m).^2);
        err_mean(m_c,c_c)=abs(m-mew);
        err_var(m_c,c_c)=abs(v-sigma^2);
        if m_c==2
            fprintf('%s, znum=%i: mean matched exactly [T1], this should be zero: %2.8e \n',nmlist{m_c},znum,err_mean(m_c,c_c))
            fprintf('%s, znum=%i: variance matched exactly [T1], this should be zero: %2.8e \n',nmlist{m_c},znum,err_var(m_c,c_c))
        elseif m_c==3
            fprintf('%s, znum=%i: mean matched to the entropy solver tolerance, this should be below %g: %2.8e \n',nmlist{m_c},znum,entropytol,err_mean(m_c,c_c))
            fprintf('%s, znum=%i: variance matched to the entropy solver tolerance, this should be below %g: %2.8e \n',nmlist{m_c},znum,entropytol,err_var(m_c,c_c))
        elseif m_c==4
            % Tauchen-Hussey places Gauss-Hermite nodes, and n-node Gauss-Hermite integrates
            % polynomials of degree up to 2n-1 exactly. At rho=0 the transition density is normal,
            % so the mean (degree 1) and variance (degree 2) come out EXACT for any znum>=2.
            fprintf('%s, znum=%i: mean exact by Gauss-Hermite exactness [T1], this should be zero: %2.8e \n',nmlist{m_c},znum,err_mean(m_c,c_c))
            fprintf('%s, znum=%i: variance exact by Gauss-Hermite exactness [T1], this should be zero: %2.8e \n',nmlist{m_c},znum,err_var(m_c,c_c))
        else
            fprintf('%s, znum=%i: mean error [T2] %2.3e, variance error [T2] %2.3e \n',nmlist{m_c},znum,err_mean(m_c,c_c),err_var(m_c,c_c))
        end

        output.method(m_c).name=nmlist{m_c};
        output.method(m_c).sweep(c_c).znum=znum;
        output.method(m_c).sweep(c_c).z_grid=z_grid;
        output.method(m_c).sweep(c_c).pi_z=pi_z;
    end
end
output.err_mean=err_mean;
output.err_var=err_var;
output.runtime=runtime;
output.nmlist=nmlist;
output.znums=znums;

%% Convergence, for the one method here that genuinely approximates
% Only AR1_Tauchen approximates the moments at rho=0. Rouwenhorst matches them algebraically,
% Farmer-Toda targets them, and Tauchen-Hussey gets them from Gauss-Hermite exactness - so asking
% any of those three to "converge" is asking something already exact to get more exact.
% And AR1_Tauchen converges only to its truncation floor at fixed Tauchen_q, not to zero.
fprintf('AR1_Tauchen variance error falls over the first half of the sweep [T2], this should be one: %i \n',err_var(1,4)<err_var(1,1))

%% Figure
figure(figure_c)
subplot(1,2,1)
plot(znums,max(err_var(1,:),10^(-18)),'o-',znums,max(err_var(2,:),10^(-18)),'s-',znums,max(err_var(3,:),10^(-18)),'^-',znums,max(err_var(4,:),10^(-18)),'v-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('AR(1) methods at rho=0: |variance error|'); xlabel('znum'); legend(nmlist,'Location','best','Interpreter','none')
subplot(1,2,2)
plot(znums,runtime(1,:),'o-',znums,runtime(2,:),'s-',znums,runtime(3,:),'^-',znums,runtime(4,:),'v-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('runtime (s)'); xlabel('znum'); legend(nmlist,'Location','best','Interpreter','none')

end
