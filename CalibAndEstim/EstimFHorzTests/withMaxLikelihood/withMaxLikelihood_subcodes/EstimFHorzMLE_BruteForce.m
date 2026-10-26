function output=EstimFHorzMLE_BruteForce()
% doPart(1): LogLikelihood_FHorz against an enumeration of every latent path, on tiny models (N_j=4) with a HAND-BUILT Policy.
%
% Nothing here comes from the toolkit except the LogLikelihood_FHorz call being tested: the Policy is written down by a formula
% (so no value function iteration), the model values of the observed variables, the next-period lotteries of the grid
% interpolation layer and the emission densities are all computed by hand in plain loops, and the likelihood is computed by
% EstimFHorzMLE_Enumerate, which never merges two paths that reach the same state (the forward algorithm does nothing else).
% So each check is an identity between two routes to the same number.
%
% Configurations (the Data ages have gaps, late entrants, an individual with one age, and missing observations in all of them):
%   (a) one endogenous state, no z, no e, no d; normal additive measurement error
%   (b) z; normal additive, normal log (multiplicative), and z observed exactly
%   (c) e, grid interpolation layer; Student-t additive, and the endogenous state observed via gridlottery
%   (d) d, z and e, grid interpolation layer; d misclassified (uniform, K=3, values given), Student-t, normal with an age-dependent sigma
%   (e) two endogenous states and z, grid interpolation layer; a1 via gridlottery, z misclassified (K from the model values), normal age-dependent sigma
%   (f) as (d), with estimoptions.initialcondition='conditional', the Data.agej0 form of the ages, and two individuals per batch

N_j=4;
Params.r=0.05;
Params.w=1;
Params.kappa_j=[1,1.2,1.4,0.8];
Params.sig_c=0.3;
Params.sig_c_j=[0.2,0.3,0.4,0.5];
Params.sig_y=0.15;
Params.nu=5;
Params.pmis=0.1;
Params.pmisz=0.2;

configs={'a','b','c','d','e','f'};
for cc=1:length(configs)
    config=configs{cc};
    clear FnsToEvaluate MeasurementError Data
    simoptions=struct();
    estimoptions=struct();
    n_d=0; d_grid=[];
    n_z=0; z_grid=[]; pi_z=[];
    n_e=0; e_grid=[]; pi_e=[];
    GI=0; ngi=0;
    conditional=0;
    switch config
        case 'a'
            n_a=6;
        case 'b'
            n_a=6;
            n_z=3; z_grid=[0.8;1;1.25]; pi_z=[0.7,0.2,0.1; 0.2,0.6,0.2; 0.1,0.3,0.6];
        case 'c'
            n_a=6;
            n_e=2; e_grid=[0.9;1.1]; pi_e=[0.4;0.6];
            GI=1; ngi=3;
        case {'d','f'}
            n_a=5;
            n_z=2; z_grid=[0.9;1.1]; pi_z=[0.8,0.2; 0.3,0.7];
            n_e=2; e_grid=[0.95;1.05]; pi_e=[0.5;0.5];
            n_d=3; d_grid=[0;0.5;1];
            GI=1; ngi=2;
            if strcmp(config,'f')
                conditional=1;
                estimoptions.initialcondition='conditional';
                estimoptions.nindividualsbatch=2;
            end
        case 'e'
            n_a=[4,3];
            n_z=3; z_grid=[0.8;1;1.25]; pi_z=[0.6,0.3,0.1; 0.25,0.5,0.25; 0.1,0.3,0.6];
            GI=1; ngi=3;
    end
    l_a=length(n_a);
    a1_grid=4*linspace(0,1,n_a(1))'.^2;
    if l_a==1
        a_grid=a1_grid;
        a2_grid=0;
    else
        a2_grid=[0;0.5;1];
        a_grid=[a1_grid;a2_grid];
    end
    N_a=prod(n_a);
    N_z=prod(n_z); N_zi=max(N_z,1);
    N_e=prod(n_e); N_ei=max(N_e,1);
    N_s=N_a*N_zi*N_ei;
    l_d=(n_d>0);
    if N_e>0
        simoptions.n_e=n_e; simoptions.e_grid=e_grid; simoptions.pi_e=pi_e;
    end
    if GI==1
        simoptions.gridinterplayer=1; simoptions.ngridinterp=ngi;
    end
    zg=ones(N_zi,1); if N_z>0, zg=z_grid; end
    eg=ones(N_ei,1); if N_e>0, eg=e_grid; end

    %% The hand-built Policy, the next-period lotteries, and the model values
    nrows=l_d+l_a+2*GI;
    PolicyH=zeros(nrows,N_a,N_zi,N_ei,N_j);
    AprimeIdx=ones(N_s,1+GI,N_j-1);
    AprimeProb=zeros(N_s,1+GI,N_j-1);
    V_cons=zeros(N_s,N_j); V_income=zeros(N_s,N_j); V_z=zeros(N_s,N_j); V_a1=zeros(N_s,N_j); V_hours=zeros(N_s,N_j); V_earn=zeros(N_s,N_j);
    for jj=1:N_j
        for ie=1:N_ei
            for iz=1:N_zi
                for ia=1:N_a
                    if l_a==1
                        ia1=ia; ia2=1;
                    else
                        [ia1,ia2]=ind2sub(n_a,ia);
                    end
                    s=ia+N_a*(iz-1)+N_a*N_zi*(ie-1);
                    rr=0;
                    if l_d==1
                        id=mod(ia1+iz+ie+jj,n_d)+1;
                        dval=d_grid(id);
                        PolicyH(1,ia,iz,ie,jj)=id;
                        rr=1;
                    else
                        dval=1;
                    end
                    if GI==1
                        a1low=mod(ia1+2*iz+3*ie+jj,n_a(1)-1)+1;
                        L2=mod(ia1+iz+2*ie+jj,ngi+1)+1;
                        wup=(L2-1)/(ngi+1);
                        a1pval=a1_grid(a1low)+wup*(a1_grid(a1low+1)-a1_grid(a1low));
                        PolicyH(rr+1,ia,iz,ie,jj)=a1low;
                        PolicyH(end-1,ia,iz,ie,jj)=L2;
                        PolicyH(end,ia,iz,ie,jj)=2; % L2flag: the usual (no forcing)
                    else
                        a1low=mod(ia1+2*iz+3*ie+jj,n_a(1))+1;
                        wup=0;
                        a1pval=a1_grid(a1low);
                        PolicyH(rr+1,ia,iz,ie,jj)=a1low;
                    end
                    if l_a==2
                        a2p=mod(ia2+iz+jj,n_a(2))+1;
                        a2pval=a2_grid(a2p);
                        PolicyH(rr+2,ia,iz,ie,jj)=a2p;
                    else
                        a2p=1;
                        a2pval=0;
                    end
                    if jj<N_j
                        AprimeIdx(s,1,jj)=a1low+n_a(1)*(a2p-1);
                        AprimeProb(s,1,jj)=1-wup;
                        if GI==1
                            AprimeIdx(s,2,jj)=a1low+1+n_a(1)*(a2p-1);
                            AprimeProb(s,2,jj)=wup;
                        end
                    end
                    aval=a1_grid(ia1)+(l_a==2)*a2_grid(ia2);
                    earn=Params.w*Params.kappa_j(jj)*zg(iz)*eg(ie)*dval;
                    V_cons(s,jj)=(1+Params.r)*aval+earn-a1pval-a2pval;
                    V_income(s,jj)=earn;
                    V_earn(s,jj)=earn;
                    V_z(s,jj)=zg(iz);
                    V_a1(s,jj)=a1_grid(ia1);
                    V_hours(s,jj)=dval;
                end
            end
        end
    end
    temp=[n_a,n_z(N_z>0),n_e(N_e>0),N_j];
    PolicyH=gpuArray(reshape(PolicyH,[nrows,temp]));

    %% The distribution at age 1: mass on the two lowest a1 points (lowest a2), every z and e, unequal weights
    mu1=zeros(N_s,1);
    cnt=0;
    for ie=1:N_ei
        for iz=1:N_zi
            for ia1=1:2
                cnt=cnt+1;
                mu1(ia1+N_a*(iz-1)+N_a*N_zi*(ie-1))=cnt;
            end
        end
    end
    mu1=mu1/sum(mu1);
    temp=[n_a,n_z(N_z>0),n_e(N_e>0)];
    if isscalar(temp)
        temp=[temp,1];
    end
    jequaloneDist=gpuArray(reshape(mu1,temp));

    pi_z_J=ones(1,1,N_j-1); if N_z>0, pi_z_J=repmat(pi_z,1,1,N_j-1); end
    pi_e_J=ones(1,N_j); if N_e>0, pi_e_J=repmat(pi_e,1,N_j); end

    %% The data
    if strcmp(config,'f')
        agej0=[1;2;1;2;2;1];
        T=3;
        Age=agej0+(0:1:T-1);
    else
        Age=[1,2,3,4; 2,3,4,NaN; 1,NaN,3,4; 3,4,NaN,NaN; 4,NaN,NaN,NaN; 1,2,NaN,4];
        T=4;
    end
    N=size(Age,1);
    switch config
        case 'a'
            obsnames={'cons'};
            FnsToEvaluate.cons=@(aprime,a,r,w,kappa_j) (1+r)*a+w*kappa_j-aprime;
            MeasurementError.cons.distribution='normal'; MeasurementError.cons.sigma='sig_c';
        case 'b'
            obsnames={'cons','income','zobs'};
            FnsToEvaluate.cons=@(aprime,a,z,r,w,kappa_j) (1+r)*a+w*kappa_j*z-aprime;
            FnsToEvaluate.income=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
            FnsToEvaluate.zobs=@(aprime,a,z) z;
            MeasurementError.cons.sigma='sig_c'; % distribution and type left at their defaults (normal, additive)
            MeasurementError.income.distribution='normal'; MeasurementError.income.type='log'; MeasurementError.income.sigma='sig_y';
            MeasurementError.zobs.distribution='exact';
        case 'c'
            obsnames={'cons','assets'};
            FnsToEvaluate.cons=@(aprime,a,e,r,w,kappa_j) (1+r)*a+w*kappa_j*e-aprime;
            FnsToEvaluate.assets=@(aprime,a,e) a;
            MeasurementError.cons.distribution='t'; MeasurementError.cons.type='additive'; MeasurementError.cons.sigma='sig_c'; MeasurementError.cons.dof='nu';
            MeasurementError.assets.distribution='gridlottery';
        case {'d','f'}
            obsnames={'hours','earnings','cons'};
            FnsToEvaluate.hours=@(d,aprime,a,z,e) d;
            FnsToEvaluate.earnings=@(d,aprime,a,z,e,w,kappa_j) w*kappa_j*z*e*d;
            FnsToEvaluate.cons=@(d,aprime,a,z,e,r,w,kappa_j) (1+r)*a+w*kappa_j*z*e*d-aprime;
            MeasurementError.hours.distribution='misclassified'; MeasurementError.hours.misclassification='uniform'; MeasurementError.hours.misclassprob='pmis'; MeasurementError.hours.values=d_grid;
            MeasurementError.earnings.distribution='t'; MeasurementError.earnings.sigma='sig_y'; MeasurementError.earnings.dof='nu';
            MeasurementError.cons.distribution='normal'; MeasurementError.cons.sigma='sig_c_j';
        case 'e'
            obsnames={'a1obs','cons','zobs'};
            FnsToEvaluate.a1obs=@(aprime1,aprime2,a1,a2,z) a1;
            FnsToEvaluate.cons=@(aprime1,aprime2,a1,a2,z,r,w,kappa_j) (1+r)*(a1+a2)+w*kappa_j*z-aprime1-aprime2;
            FnsToEvaluate.zobs=@(aprime1,aprime2,a1,a2,z) z;
            MeasurementError.a1obs.distribution='gridlottery';
            MeasurementError.cons.distribution='normal'; MeasurementError.cons.sigma='sig_c_j';
            MeasurementError.zobs.distribution='misclassified'; MeasurementError.zobs.misclassprob='pmisz'; % no .values: K from the model values (=3)
    end
    nObs=length(obsnames);
    ObsByAgeH=NaN(N,N_j,nObs);
    for vv=1:nObs
        Data.(obsnames{vv})=NaN(N,T);
        switch obsnames{vv}
            case 'cons'
                Vv=V_cons;
            case 'income'
                Vv=V_income;
            case 'earnings'
                Vv=V_earn;
            case 'zobs'
                Vv=V_z;
            case {'assets','a1obs'}
                Vv=V_a1;
            case 'hours'
                Vv=V_hours;
        end
        distv='normal';
        if isfield(MeasurementError.(obsnames{vv}),'distribution')
            distv=MeasurementError.(obsnames{vv}).distribution;
        end
        for ii=1:N
            % gridlottery: the observation at each age is put next to an a1 grid point that the individual can actually be at,
            % given the gridlottery observations at earlier ages (Sind is the set of such states; the other observations of
            % these configurations have positive emission at every state, so they never shrink it)
            Sind=find(mu1>0);
            for jj=1:N_j
                tt=find(Age(ii,:)==jj);
                if strcmp(distv,'gridlottery') && ~isempty(tt)
                    tempa1=mod(Sind-1,n_a(1))+1; % a1 index of the possible states
                    tempa1=unique(tempa1(tempa1<n_a(1)));
                    kk=tempa1(mod(ii+tt,numel(tempa1))+1);
                    frac=0.2+0.6*mod(ii+tt,2);
                    Data.(obsnames{vv})(ii,tt)=a1_grid(kk)+frac*(a1_grid(kk+1)-a1_grid(kk));
                    tempa1=mod(Sind-1,n_a(1))+1;
                    Sind=Sind(tempa1==kk | tempa1==kk+1);
                end
                if jj<N_j
                    temp=[];
                    for pp=1:(1+GI)
                        ap=AprimeIdx(Sind,pp,jj);
                        ap=ap(AprimeProb(Sind,pp,jj)>0);
                        for zp=1:N_zi
                            for ep=1:N_ei
                                temp=[temp; ap+N_a*(zp-1)+N_a*N_zi*(ep-1)];
                            end
                        end
                    end
                    Sind=unique(temp);
                end
            end
            for tt=1:T
                jj=Age(ii,tt);
                if isnan(jj)
                    continue
                end
                switch distv
                    case {'normal','t'}
                        if isfield(MeasurementError.(obsnames{vv}),'type') && strcmp(MeasurementError.(obsnames{vv}).type,'log')
                            y=exp(mean(log(Vv(:,jj))))*exp(0.1*sin(3*ii+2*tt+vv));
                        else
                            y=mean(Vv(:,jj))+0.5*std(Vv(:,jj))*sin(3*ii+2*tt+vv);
                        end
                    case 'exact'
                        y=zg(mod(ii+tt,N_zi)+1);
                    case 'misclassified'
                        if strcmp(obsnames{vv},'hours')
                            y=d_grid(mod(2*ii+tt,n_d)+1);
                        else
                            y=zg(mod(ii+2*tt,N_zi)+1);
                        end
                    case 'gridlottery'
                        y=Data.(obsnames{vv})(ii,tt); % already set above
                end
                Data.(obsnames{vv})(ii,tt)=y;
            end
        end
        % missing observations
        Data.(obsnames{vv})(1,2)=NaN;
        Data.(obsnames{vv})(6,1)=NaN; % individual 6 enters at age 1 with nothing observed there
    end
    if strcmp(config,'f')
        Data.agej0=agej0;
    else
        Data.age=Age;
    end
    for vv=1:nObs
        for ii=1:N
            for tt=1:T
                if ~isnan(Age(ii,tt))
                    ObsByAgeH(ii,Age(ii,tt),vv)=Data.(obsnames{vv})(ii,tt);
                end
            end
        end
    end
    j0=min(Age,[],2);
    jlast=max(Age,[],2);

    %% The log emissions by hand
    logEmis=zeros(N_s,N,N_j);
    for vv=1:nObs
        ME=MeasurementError.(obsnames{vv});
        switch obsnames{vv}
            case 'cons'
                Vv=V_cons;
            case 'income'
                Vv=V_income;
            case 'earnings'
                Vv=V_earn;
            case 'zobs'
                Vv=V_z;
            case {'assets','a1obs'}
                Vv=V_a1;
            case 'hours'
                Vv=V_hours;
        end
        if ~isfield(ME,'distribution')
            ME.distribution='normal';
        end
        for ii=1:N
            for jj=1:N_j
                y=ObsByAgeH(ii,jj,vv);
                if isnan(y)
                    continue
                end
                for s=1:N_s
                    m=Vv(s,jj);
                    switch ME.distribution
                        case 'normal'
                            sg=Params.(ME.sigma); if numel(sg)>1, sg=sg(jj); end
                            if isfield(ME,'type') && strcmp(ME.type,'log')
                                f=exp(-0.5*((log(y)-log(m))/sg)^2)/(sg*y*sqrt(2*pi));
                            else
                                f=exp(-0.5*((y-m)/sg)^2)/(sg*sqrt(2*pi));
                            end
                        case 't'
                            sg=Params.(ME.sigma); if numel(sg)>1, sg=sg(jj); end
                            nu=Params.(ME.dof);
                            f=gamma((nu+1)/2)/(gamma(nu/2)*sqrt(nu*pi)*sg)*(1+((y-m)/sg)^2/nu)^(-(nu+1)/2);
                        case 'exact'
                            f=double(y==m);
                        case 'misclassified'
                            p=Params.(ME.misclassprob);
                            K=3; % hours: d_grid has 3 values; zobs: z_grid has 3 values
                            if y==m
                                f=1-p;
                            else
                                f=p/(K-1);
                            end
                        case 'gridlottery'
                            kk=find(a1_grid<=y,1,'last');
                            wgt=(y-a1_grid(kk))/(a1_grid(kk+1)-a1_grid(kk));
                            if m==a1_grid(kk)
                                f=1-wgt;
                            elseif m==a1_grid(kk+1)
                                f=wgt;
                            else
                                f=0;
                            end
                    end
                    logEmis(s,ii,jj)=logEmis(s,ii,jj)+log(f);
                end
            end
        end
    end

    %% The two routes
    [LogLijBF,LogLiBF,npaths]=EstimFHorzMLE_Enumerate(mu1,AprimeIdx,AprimeProb,pi_z_J,pi_e_J,logEmis,N_a,N_zi,N_ei,N_j,j0,jlast,conditional);
    [LogL,LogLi,LogLij]=LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,PolicyH,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
    LogLij=gather(LogLij); LogLi=gather(LogLi); LogL=gather(LogL);

    fin=~isnan(LogLijBF) & ~isnan(LogLij);
    fprintf('BruteForce (%s): %i paths enumerated, LogL=%.6f, smallest per-age term %.4f (informational) \n',config,npaths,LogL,min(LogLijBF(:)))
    fprintf('BruteForce (%s) LogLij vs the enumeration of every path, this should be zero: %.3e, NaN-pattern mismatches, this should be zero: %i \n',config,max([0;abs(LogLij(fin)-LogLijBF(fin))]),sum(isnan(LogLij(:))~=isnan(LogLijBF(:))))
    fprintf('BruteForce (%s) LogLi vs the enumeration of every path, this should be zero: %.3e \n',config,CoreMaxAbs(LogLi-LogLiBF))
    fprintf('BruteForce (%s) LogL vs sum(LogLi), this should be zero: %.3e \n',config,CoreMaxAbs(LogL-sum(LogLiBF)))
    output.(['LogLi_',config])=LogLi;
    output.(['LogLiBF_',config])=LogLiBF;
end

end
