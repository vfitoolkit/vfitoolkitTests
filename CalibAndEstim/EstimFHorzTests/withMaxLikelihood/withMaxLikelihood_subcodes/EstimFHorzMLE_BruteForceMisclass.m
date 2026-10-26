function output=EstimFHorzMLE_BruteForceMisclass()
% doPart(9): LogLikelihood_FHorz with the misclassification schemes 'matrix', 'ordered' and 'frequency' (model and data shares),
% against an enumeration of every latent path (EstimFHorzMLE_Enumerate) on a tiny model (N_j=4) with a HAND-BUILT Policy (as doPart(1)).
% The model has a discrete decision d (hours, three values), z (three values) and the grid interpolation layer; the observed
% variables are hours and z (both misclassified) and consumption (normal measurement error).
% The misclassification matrices E(k,m)=P(report k | true m) are built here by hand; the model's shares of hours at each age (for
% 'frequency' with model shares) from a hand forward iteration of the population distribution (no observations).
%   (j) hours 'matrix' (a different matrix at each age, K-by-K-by-N_j), z 'ordered'
%   (k) hours 'frequency' with the model's shares, z 'frequency' with the data's shares; the conditional initial condition
%   (l) hours 'uniform' and the same as a 'matrix': the two give the same LogLi, and both equal the enumeration
%   (m) refusals

N_j=4;
Params.r=0.05; Params.w=1; Params.kappa_j=[1,1.2,1.4,0.8];
Params.sig_c=0.3;
Params.sig_ord=0.15; % 'ordered' s.d. (z values are 0.8, 1, 1.25, so the midpoints are 0.1 and 0.125 away)
Params.p_h=0.2; Params.p_z=0.3; Params.pmis=0.15;
Ehours=zeros(3,3,N_j);
for jj=1:N_j
    temp=[0.8,0.1,0.05; 0.15,0.8,0.15; 0.05,0.1,0.8]+0.02*jj*[0,1,0;0,-1,0;0,0,0];
    Ehours(:,:,jj)=temp./sum(temp,1);
end
Params.Ehours=Ehours;
K=3;
Euniform=(1-Params.pmis)*eye(K)+(Params.pmis/(K-1))*(1-eye(K));
Params.Euniform=Euniform;

n_d=3; d_grid=[0;0.5;1];
n_a=4; a_grid=4*linspace(0,1,n_a)'.^2;
n_z=3; z_grid=[0.8;1;1.25]; pi_z=[0.7,0.2,0.1; 0.2,0.6,0.2; 0.1,0.3,0.6];
ngi=2;
N_a=n_a; N_z=n_z; N_s=N_a*N_z;
simoptions.gridinterplayer=1; simoptions.ngridinterp=ngi;

% The hand-built Policy, the lotteries, the model values (hours index, z index, consumption)
PolicyH=zeros(4,N_a,N_z,N_j);
AprimeIdx=ones(N_s,2,N_j-1); AprimeProb=zeros(N_s,2,N_j-1);
Hidx=zeros(N_s,N_j); Zidx=zeros(N_s,N_j); Vcons=zeros(N_s,N_j);
for jj=1:N_j
    for iz=1:N_z
        for ia=1:N_a
            s=ia+N_a*(iz-1);
            id=mod(ia+iz+jj,n_d)+1;
            a1low=mod(ia+2*iz+jj,n_a-1)+1;
            L2=mod(ia+iz+jj,ngi+1)+1;
            wup=(L2-1)/(ngi+1);
            PolicyH(:,ia,iz,jj)=[id;a1low;L2;2];
            if jj<N_j
                AprimeIdx(s,:,jj)=[a1low,a1low+1]; AprimeProb(s,:,jj)=[1-wup,wup];
            end
            Hidx(s,jj)=id; Zidx(s,jj)=iz;
            Vcons(s,jj)=(1+Params.r)*a_grid(ia)+Params.w*Params.kappa_j(jj)*z_grid(iz)*d_grid(id)-(a_grid(a1low)+wup*(a_grid(a1low+1)-a_grid(a1low)));
        end
    end
end
PolicyH=gpuArray(PolicyH);
mu=zeros(N_a,N_z); mu(1:2,:)=reshape(1:6,[2,3]); mu=mu/sum(mu(:)); mu1=mu(:);
jequaloneDist=gpuArray(mu);
pi_z_J=repmat(pi_z,1,1,N_j-1);
% The population distribution at each age (hand forward iteration, no observations), and the model's shares of hours
pop=zeros(N_s,N_j); pop(:,1)=mu1;
for jj=1:N_j-1
    for s=1:N_s
        iz=floor((s-1)/N_a)+1;
        for pp=1:2
            for zp=1:N_z
                sp=AprimeIdx(s,pp,jj)+N_a*(zp-1);
                pop(sp,jj+1)=pop(sp,jj+1)+pop(s,jj)*AprimeProb(s,pp,jj)*pi_z(iz,zp);
            end
        end
    end
end
qh=zeros(K,N_j);
for jj=1:N_j
    for kk=1:K
        qh(kk,jj)=sum(pop(Hidx(:,jj)==kk,jj));
    end
end

FnsToEvaluate.hours=@(d,aprime,a,z) d;
FnsToEvaluate.zobs=@(d,aprime,a,z) z;
FnsToEvaluate.cons=@(d,aprime,a,z,r,w,kappa_j) (1+r)*a+w*kappa_j*z*d-aprime;
Age=[1,2,3,4; 2,3,4,NaN; 1,NaN,3,4; 3,4,NaN,NaN; 4,NaN,NaN,NaN; 1,2,NaN,4];
N=size(Age,1); T=4;
j0=min(Age,[],2); jlast=max(Age,[],2);

configs={'j','k','l'};
for cc=1:length(configs)
    config=configs{cc};
    clear MeasurementError
    estimoptions=struct();
    conditional=0;
    MeasurementError.cons.distribution='normal'; MeasurementError.cons.sigma='sig_c';
    MeasurementError.hours.distribution='misclassified';
    MeasurementError.zobs.distribution='misclassified';
    switch config
        case 'j'
            MeasurementError.hours.misclassification='matrix'; MeasurementError.hours.misclassmatrix='Ehours'; MeasurementError.hours.values=d_grid;
            MeasurementError.zobs.misclassification='ordered'; MeasurementError.zobs.misclasssigma='sig_ord';
        case 'k'
            MeasurementError.hours.misclassification='frequency'; MeasurementError.hours.misclassprob='p_h'; MeasurementError.hours.values=d_grid; % frequencyshares left at the default, 'model'
            MeasurementError.zobs.misclassification='frequency'; MeasurementError.zobs.misclassprob='p_z'; MeasurementError.zobs.frequencyshares='data';
            conditional=1; estimoptions.initialcondition='conditional';
        case 'l'
            MeasurementError.hours.misclassification='uniform'; MeasurementError.hours.misclassprob='pmis'; MeasurementError.hours.values=d_grid;
            MeasurementError.zobs.misclassification='matrix'; MeasurementError.zobs.misclassmatrix='Euniform';
    end
    % The data (hours only at values the model has at that age, so that 'frequency' with model shares has a positive likelihood)
    clear Data
    Data.age=Age;
    Data.hours=NaN(N,T); Data.zobs=NaN(N,T); Data.cons=NaN(N,T);
    for ii=1:N
        for tt=1:T
            jj=Age(ii,tt);
            if isnan(jj), continue, end
            posk=find(qh(:,jj)>0);
            Data.hours(ii,tt)=d_grid(posk(mod(2*ii+tt,numel(posk))+1));
            Data.zobs(ii,tt)=z_grid(mod(ii+2*tt,N_z)+1);
            Data.cons(ii,tt)=mean(Vcons(:,jj))+0.5*std(Vcons(:,jj))*sin(3*ii+2*tt);
        end
    end
    Data.hours(1,2)=NaN; Data.zobs(6,1)=NaN; Data.cons(6,1)=NaN; Data.hours(6,1)=NaN;
    % The misclassification matrices by hand, at each age
    Eh=zeros(K,K,N_j); Ez=zeros(K,K,N_j);
    for jj=1:N_j
        switch config
            case 'j'
                Eh(:,:,jj)=Ehours(:,:,jj);
                c=[-Inf;(z_grid(1:end-1)+z_grid(2:end))/2;Inf];
                for kk=1:K
                    for mm=1:K
                        Ez(kk,mm,jj)=0.5*erfc(-(c(kk+1)-z_grid(mm))/(Params.sig_ord*sqrt(2)))-0.5*erfc(-(c(kk)-z_grid(mm))/(Params.sig_ord*sqrt(2)));
                    end
                end
            case 'k'
                for kk=1:K
                    for mm=1:K
                        Eh(kk,mm,jj)=(1-Params.p_h)*(kk==mm)+Params.p_h*qh(kk,jj);
                    end
                end
                yz=[];
                for ii=1:N
                    tt=find(Age(ii,:)==jj);
                    if ~isempty(tt) && ~isnan(Data.zobs(ii,tt))
                        yz=[yz; find(z_grid==Data.zobs(ii,tt))];
                    end
                end
                qz=zeros(K,1);
                for kk=1:K
                    qz(kk)=sum(yz==kk)/max(numel(yz),1);
                end
                for kk=1:K
                    for mm=1:K
                        Ez(kk,mm,jj)=(1-Params.p_z)*(kk==mm)+Params.p_z*qz(kk);
                    end
                end
            case 'l'
                for kk=1:K
                    for mm=1:K
                        Eh(kk,mm,jj)=(kk==mm)*(1-Params.pmis)+(kk~=mm)*Params.pmis/(K-1);
                    end
                end
                Ez(:,:,jj)=Euniform;
        end
    end
    % The log emissions by hand
    logEmis=zeros(N_s,N,N_j);
    for ii=1:N
        for tt=1:T
            jj=Age(ii,tt);
            if isnan(jj), continue, end
            for s=1:N_s
                le=0;
                if ~isnan(Data.hours(ii,tt))
                    le=le+log(Eh(find(d_grid==Data.hours(ii,tt)),Hidx(s,jj),jj));
                end
                if ~isnan(Data.zobs(ii,tt))
                    le=le+log(Ez(find(z_grid==Data.zobs(ii,tt)),Zidx(s,jj),jj));
                end
                if ~isnan(Data.cons(ii,tt))
                    le=le+log(exp(-0.5*((Data.cons(ii,tt)-Vcons(s,jj))/Params.sig_c)^2)/(Params.sig_c*sqrt(2*pi)));
                end
                logEmis(s,ii,jj)=le;
            end
        end
    end
    [LogLijBF,LogLiBF,npaths]=EstimFHorzMLE_Enumerate(mu1,AprimeIdx,AprimeProb,pi_z_J,ones(1,N_j),logEmis,N_a,N_z,1,N_j,j0,jlast,conditional);
    [LogL,LogLi,LogLij]=LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,PolicyH,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
    LogLij=gather(LogLij); LogLi=gather(LogLi);
    fin=~isnan(LogLijBF) & ~isnan(LogLij);
    fprintf('BruteForceMisclass (%s): %i paths enumerated, LogL=%.6f (informational) \n',config,npaths,gather(LogL))
    fprintf('BruteForceMisclass (%s) LogLij vs the enumeration of every path, this should be zero: %.3e, NaN-pattern mismatches, this should be zero: %i \n',config,max([0;abs(LogLij(fin)-LogLijBF(fin))]),sum(isnan(LogLij(:))~=isnan(LogLijBF(:))))
    fprintf('BruteForceMisclass (%s) LogLi vs the enumeration of every path, this should be zero: %.3e \n',config,CoreMaxAbs(LogLi-LogLiBF))
    if strcmp(config,'l')
        ME2=MeasurementError; ME2.hours.misclassification='matrix'; ME2.hours=rmfield(ME2.hours,'misclassprob'); ME2.hours.misclassmatrix='Euniform';
        [~,LogLi2]=LogLikelihood_FHorz(Data,ME2,jequaloneDist,PolicyH,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
        fprintf('BruteForceMisclass (l) hours uniform vs the same as a matrix, this should be zero: %.3e \n',CoreMaxAbs(gather(LogLi2)-LogLi))
    end
end

%% (m) refusals
errcases={'unknown misclassification','matrix of the wrong size','matrix column not summing to one','matrix with a negative entry','ordered with values not ascending','frequencyshares neither model nor data','ordered without misclasssigma'};
errmsgs={'is not one of uniform, matrix, frequency, ordered','must be K-by-K or K-by-K-by-N_j','each column must sum to one','each column must sum to one','must be in ascending order','must be ''model'' or ''data''','misclasssigma is required'};
for ee=1:length(errcases)
    ME=MeasurementError; P=Params;
    ME.hours.misclassification='matrix'; ME.hours=rmfield(ME.hours,'misclassprob'); ME.hours.misclassmatrix='Ehours';
    switch ee
        case 1
            ME.hours.misclassification='poisson';
        case 2
            P.Ehours=eye(2);
        case 3
            P.Ehours=Ehours; P.Ehours(1,1,2)=P.Ehours(1,1,2)+0.1;
        case 4
            P.Ehours=[1.2,0,0; -0.2,1,0; 0,0,1];
        case 5
            ME.zobs=struct('distribution','misclassified','misclassification','ordered','misclasssigma','sig_ord','values',flipud(z_grid));
        case 6
            ME.zobs=struct('distribution','misclassified','misclassification','frequency','misclassprob','p_z','frequencyshares','population');
        case 7
            ME.zobs=struct('distribution','misclassified','misclassification','ordered');
    end
    errok=0;
    try
        LogLikelihood_FHorz(Data,ME,jequaloneDist,PolicyH,FnsToEvaluate,P,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions);
    catch err
        errok=contains(err.message,errmsgs{ee});
        if ~errok
            fprintf('   (the error was: %s) \n',err.message);
        end
    end
    fprintf('BruteForceMisclass (m) %s errors with the intended message, this should be zero: %i \n',errcases{ee},1-errok);
end

output.LogLi=LogLi;
end
