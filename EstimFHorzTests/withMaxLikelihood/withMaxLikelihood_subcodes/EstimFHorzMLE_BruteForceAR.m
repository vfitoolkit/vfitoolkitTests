function output=EstimFHorzMLE_BruteForceAR()
% doPart(6): LogLikelihood_FHorz with AR(1) measurement error components, against an enumeration of every latent path, on tiny models
% (N_j=4) with a HAND-BUILT Policy (as doPart(1)).
%
% The AR(1) component u of the measurement error is a Markov state independent of the model, so for the enumeration it is folded into
% the exogenous Markov state: the joint state (z,u) has transition kron(pi_u,pi_z) (z fastest), and the enumerator
% (EstimFHorzMLE_Enumerate, which never merges paths) runs on states (a,(z,u)). At the first age u is distributed as the stationary
% distribution of the discretized AR(1), computed here as the eigenvector of its transition matrix (the toolkit uses the binomial
% formula). The emission of a variable with an AR(1) component at state (a,z,u) is that of y-(model value)-u under the iid part
% (y/(model value*exp(u)) for type='log').
%
% Configurations (the Data ages have gaps, late entrants, an individual with one age, and missing observations):
%   (h) z, grid interpolation layer; cons with an additive AR(1) component, earnings with iid normal error
%   (i) no z, grid interpolation layer; cons with an additive AR(1) component and income with a multiplicative (log) one, so two AR(1)
%       states; the conditional initial condition

N_j=4;
Params.r=0.05;
Params.w=1;
Params.kappa_j=[1,1.2,1.4,0.8];
Params.sig_c=0.3;
Params.sig_y=0.2;
Params.rho_c=0.6;
Params.sigAR_c=0.25;
Params.rho_y=0.8;
Params.sigAR_y=0.1;
nAR=3; % points in each AR(1) discretization

configs={'h','i'};
for cc=1:length(configs)
    config=configs{cc};
    clear FnsToEvaluate MeasurementError Data
    simoptions=struct();
    estimoptions=struct();
    estimoptions.nARpoints=nAR;
    conditional=0;
    n_d=0; d_grid=[];
    GI=1; ngi=2;
    simoptions.gridinterplayer=1; simoptions.ngridinterp=ngi;
    switch config
        case 'h'
            n_a=4;
            n_z=2; z_grid=[0.9;1.1]; pi_z=[0.8,0.2; 0.3,0.7];
            ARvars={'cons'};
        case 'i'
            n_a=3;
            n_z=0; z_grid=[]; pi_z=[];
            ARvars={'cons','income'};
            conditional=1;
            estimoptions.initialcondition='conditional';
    end
    a_grid=4*linspace(0,1,n_a)'.^2;
    N_a=n_a;
    N_z=prod(n_z); N_zi=max(N_z,1);
    zg=ones(N_zi,1); if N_z>0, zg=z_grid; end
    pi_zi=1; if N_z>0, pi_zi=pi_z; end

    % The discretized AR(1) components, and their joint chain (the first AR variable fastest)
    pi_u=1; stat_u=1; uvals=struct();
    nARvars=length(ARvars);
    N_u=nAR^nARvars;
    for kk=1:nARvars
        if strcmp(ARvars{kk},'cons')
            [ug,pu]=discretizeAR1_Rouwenhorst(0,Params.rho_c,Params.sigAR_c,nAR);
        else
            [ug,pu]=discretizeAR1_Rouwenhorst(0,Params.rho_y,Params.sigAR_y,nAR);
        end
        ug=gather(ug); pu=gather(pu);
        [V,Dg]=eig(pu');
        [~,kmax]=max(real(diag(Dg)));
        su=real(V(:,kmax)); su=su/sum(su); % stationary distribution by hand
        pi_u=kron(pu,pi_u);
        stat_u=kron(su,stat_u);
        temp=zeros(N_u,1);
        for uu=1:N_u
            idx=mod(floor((uu-1)/nAR^(kk-1)),nAR)+1; % this component's index at joint index uu
            temp(uu)=ug(idx);
        end
        uvals.(ARvars{kk})=temp;
    end
    N_zu=N_zi*N_u;
    N_sE=N_a*N_zu; % enumerator states (a,(z,u))

    %% The hand-built Policy, the next-period lotteries (on the enumerator states), and the model values (on (a,z))
    PolicyH=zeros(3,N_a,N_zi,N_j);
    AprimeIdx=ones(N_sE,2,N_j-1);
    AprimeProb=zeros(N_sE,2,N_j-1);
    V_cons=zeros(N_a,N_zi,N_j); V_earn=zeros(N_a,N_zi,N_j);
    for jj=1:N_j
        for iz=1:N_zi
            for ia=1:N_a
                a1low=mod(ia+2*iz+jj,n_a-1)+1;
                L2=mod(ia+iz+jj,ngi+1)+1;
                wup=(L2-1)/(ngi+1);
                a1pval=a_grid(a1low)+wup*(a_grid(a1low+1)-a_grid(a1low));
                PolicyH(1,ia,iz,jj)=a1low;
                PolicyH(2,ia,iz,jj)=L2;
                PolicyH(3,ia,iz,jj)=2; % L2flag: the usual (no forcing)
                earn=Params.w*Params.kappa_j(jj)*zg(iz);
                V_cons(ia,iz,jj)=(1+Params.r)*a_grid(ia)+earn-a1pval;
                V_earn(ia,iz,jj)=earn;
                if jj<N_j
                    for uu=1:N_u
                        s=ia+N_a*(iz-1+N_zi*(uu-1));
                        AprimeIdx(s,1,jj)=a1low; AprimeProb(s,1,jj)=1-wup;
                        AprimeIdx(s,2,jj)=a1low+1; AprimeProb(s,2,jj)=wup;
                    end
                end
            end
        end
    end
    temp=[n_a,n_z(N_z>0),N_j];
    PolicyH=gpuArray(reshape(PolicyH,[3,temp]));

    % The distribution at age 1: mass on the two lowest a points, every z, unequal weights; u at its stationary distribution
    mu_s=zeros(N_a,N_zi);
    cnt=0;
    for iz=1:N_zi
        for ia=1:2
            cnt=cnt+1;
            mu_s(ia,iz)=cnt;
        end
    end
    mu_s=mu_s/sum(mu_s(:));
    mu1=kron(stat_u,mu_s(:)); % (a,z) fastest, then u
    temp=[n_a,n_z(N_z>0)];
    if isscalar(temp)
        temp=[temp,1];
    end
    jequaloneDist=gpuArray(reshape(mu_s,temp));
    pi_zu_J=repmat(kron(pi_u,pi_zi),1,1,N_j-1);
    pi_e_J=ones(1,N_j);

    %% The data
    Age=[1,2,3,4; 2,3,4,NaN; 1,NaN,3,4; 3,4,NaN,NaN; 4,NaN,NaN,NaN; 1,2,NaN,4];
    T=4; N=size(Age,1);
    if strcmp(config,'h')
        obsnames={'cons','earnings'};
        FnsToEvaluate.cons=@(aprime,a,z,r,w,kappa_j) (1+r)*a+w*kappa_j*z-aprime;
        FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
        MeasurementError.cons.distribution='normal'; MeasurementError.cons.sigma='sig_c'; MeasurementError.cons.rho='rho_c'; MeasurementError.cons.sigmaAR='sigAR_c';
        MeasurementError.earnings.distribution='normal'; MeasurementError.earnings.sigma='sig_y';
        Vlist={V_cons,V_earn};
    else
        obsnames={'cons','income'};
        FnsToEvaluate.cons=@(aprime,a,r,w,kappa_j) (1+r)*a+w*kappa_j-aprime;
        FnsToEvaluate.income=@(aprime,a,w,kappa_j) w*kappa_j;
        MeasurementError.cons.distribution='normal'; MeasurementError.cons.sigma='sig_c'; MeasurementError.cons.rho='rho_c'; MeasurementError.cons.sigmaAR='sigAR_c';
        MeasurementError.income.distribution='normal'; MeasurementError.income.type='log'; MeasurementError.income.sigma='sig_y'; MeasurementError.income.rho='rho_y'; MeasurementError.income.sigmaAR='sigAR_y';
        Vlist={V_cons,V_earn};
    end
    nObs=length(obsnames);
    ObsByAgeH=NaN(N,N_j,nObs);
    for vv=1:nObs
        Vv=reshape(Vlist{vv},[N_a*N_zi,N_j]);
        Data.(obsnames{vv})=NaN(N,T);
        for ii=1:N
            for tt=1:T
                jj=Age(ii,tt);
                if isnan(jj)
                    continue
                end
                if strcmp(obsnames{vv},'income')
                    Data.(obsnames{vv})(ii,tt)=exp(mean(log(Vv(:,jj))))*exp(0.15*sin(3*ii+2*tt+vv));
                else
                    Data.(obsnames{vv})(ii,tt)=mean(Vv(:,jj))+0.5*std(Vv(:,jj))*sin(3*ii+2*tt+vv);
                end
            end
        end
        Data.(obsnames{vv})(1,2)=NaN;
        Data.(obsnames{vv})(6,1)=NaN;
        for ii=1:N
            for tt=1:T
                if ~isnan(Age(ii,tt))
                    ObsByAgeH(ii,Age(ii,tt),vv)=Data.(obsnames{vv})(ii,tt);
                end
            end
        end
    end
    Data.age=Age;
    j0=min(Age,[],2);
    jlast=max(Age,[],2);

    %% The log emissions by hand, on the enumerator states (a,(z,u))
    logEmis=zeros(N_sE,N,N_j);
    for vv=1:nObs
        ME=MeasurementError.(obsnames{vv});
        sg=Params.(ME.sigma);
        Vv=reshape(Vlist{vv},[N_a*N_zi,N_j]);
        for ii=1:N
            for jj=1:N_j
                y=ObsByAgeH(ii,jj,vv);
                if isnan(y)
                    continue
                end
                for uu=1:N_u
                    if isfield(ME,'rho')
                        uv=uvals.(obsnames{vv})(uu);
                    else
                        uv=0;
                    end
                    for sa=1:N_a*N_zi
                        m=Vv(sa,jj);
                        if isfield(ME,'type') && strcmp(ME.type,'log')
                            f=exp(-0.5*((log(y)-log(m)-uv)/sg)^2)/(sg*y*sqrt(2*pi));
                        else
                            f=exp(-0.5*((y-m-uv)/sg)^2)/(sg*sqrt(2*pi));
                        end
                        s=sa+N_a*N_zi*(uu-1);
                        logEmis(s,ii,jj)=logEmis(s,ii,jj)+log(f);
                    end
                end
            end
        end
    end

    %% The two routes
    [LogLijBF,LogLiBF,npaths]=EstimFHorzMLE_Enumerate(mu1,AprimeIdx,AprimeProb,pi_zu_J,pi_e_J,logEmis,N_a,N_zu,1,N_j,j0,jlast,conditional);
    [LogL,LogLi,LogLij]=LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,PolicyH,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
    LogLij=gather(LogLij); LogLi=gather(LogLi); LogL=gather(LogL);

    fin=~isnan(LogLijBF) & ~isnan(LogLij);
    fprintf('BruteForceAR (%s): %i paths enumerated, LogL=%.6f (informational) \n',config,npaths,LogL)
    fprintf('BruteForceAR (%s) LogLij vs the enumeration of every path, this should be zero: %.3e, NaN-pattern mismatches, this should be zero: %i \n',config,max([0;abs(LogLij(fin)-LogLijBF(fin))]),sum(isnan(LogLij(:))~=isnan(LogLijBF(:))))
    fprintf('BruteForceAR (%s) LogLi vs the enumeration of every path, this should be zero: %.3e \n',config,CoreMaxAbs(LogLi-LogLiBF))
    output.(['LogLi_',config])=LogLi;
end

end
