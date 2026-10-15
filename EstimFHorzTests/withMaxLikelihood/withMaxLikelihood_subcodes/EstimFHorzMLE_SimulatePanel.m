function Data=EstimFHorzMLE_SimulatePanel(Policy,n_a,n_z,N_j,a_grid,z_grid,pi_z,jequaloneDist,Params,ngridinterp,N,seed)
% A hand-coded simulation of one panel from the EstimFHorzTests model (one endogenous state, z, grid interpolation layer), at the
% parameters in Params, for the maximum likelihood tests. Plain loops on the cpu with a fixed seed, so the same panel every time.
% (Test-bank only: the toolkit's estimation never simulates.)
%
% Every individual starts at age 1 (drawn from jequaloneDist) and is followed to age N_j. Observed, with additive normal measurement error:
%   cons (s.d. Params.sig_c), at every age, except a few missing
%   earnings (s.d. Params.sig_y), at every age
%   assets (s.d. Params.sig_a), only at ages 5, 10 and 15
% The model values mirror the FnsToEvaluate of the tests: cons=(1+r)*a+(agej<Jr)*w*kappa_j*z+(agej>=Jr)*pension-aprime, with aprime
% the grid-interpolation-layer value; earnings=(agej<Jr)*w*kappa_j*z; assets=a.

rng(seed,'twister');
Policy=gather(Policy);
a_grid=gather(a_grid);
z_grid=gather(z_grid);
pi_z=gather(pi_z);
mu1=gather(jequaloneDist(:));
cmu1=cumsum(mu1);
cpi_z=cumsum(pi_z,2);

Data.agej0=ones(N,1);
Data.cons=NaN(N,N_j);
Data.earnings=NaN(N,N_j);
Data.assets=NaN(N,N_j);
for ii=1:N
    s=find(cmu1>=rand,1,'first');
    ia=mod(s-1,n_a)+1;
    iz=floor((s-1)/n_a)+1;
    for jj=1:N_j
        a1low=Policy(1,ia,iz,jj);
        L2=Policy(2,ia,iz,jj);
        L2flag=Policy(3,ia,iz,jj);
        % value of aprime (as PolicyInd2Val: the fine grid point, no flag forcing)
        if a1low<n_a
            aprimeval=a_grid(a1low)+(L2-1)/(ngridinterp+1)*(a_grid(a1low+1)-a_grid(a1low));
        else
            aprimeval=a_grid(a1low);
        end
        earn=(Params.agej(jj)<Params.Jr)*Params.w*Params.kappa_j(jj)*z_grid(iz);
        cons=(1+Params.r)*a_grid(ia)+earn+(Params.agej(jj)>=Params.Jr)*Params.pension-aprimeval;
        Data.cons(ii,jj)=cons+Params.sig_c*randn;
        Data.earnings(ii,jj)=earn+Params.sig_y*randn;
        if any(jj==[5,10,15])
            Data.assets(ii,jj)=a_grid(ia)+Params.sig_a*randn;
        end
        % next period: aprime lottery (as the agent distribution, with the flag forcing), then zprime
        if jj<N_j
            pup=(L2-1)/(ngridinterp+1);
            if L2flag==1
                pup=0;
            elseif L2flag==3
                pup=1;
            end
            if rand<pup
                ia=a1low+1;
            else
                ia=a1low;
            end
            iz=find(cpi_z(iz,:)>=rand,1,'first');
        end
    end
end
% a few missing consumption observations
missing=(rand(N,N_j)<0.05);
Data.cons(missing)=NaN;

end
