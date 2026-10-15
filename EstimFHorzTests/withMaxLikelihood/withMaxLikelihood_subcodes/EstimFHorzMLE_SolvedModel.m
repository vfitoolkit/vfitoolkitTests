function output=EstimFHorzMLE_SolvedModel(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions)
% doPart(2): LogLikelihood_FHorz on a solved model, against other toolkit commands and against itself. Two models, each runs (a)-(g):
%   variant 'z':  the EstimFHorzTests model as passed in (z, grid interpolation layer)
%   variant 'ze': the same with an iid e added (CoreFHorzTests e grid, ReturnFn_nod_z_e_nosemiz, jequaloneDist spread over e by pi_e)
% In both the age weights are replaced by NON-uniform ones (the EstimFHorzTests ones are uniform), so (a) can tell weighted from not.
%   (a) PopDist (the filter's population column, no observations) times the age weights equals StationaryDist_FHorz_Case1
%   (b) a cross-section (T=1, every individual at one age) equals the closed-form mixture
%       log( sum_s StationaryDist(s,j)/mewj(j) * prod_v f_v(y_iv | model value at s) ), with the emission densities by hand
%       (normal additive, Student-t additive, and gridlottery for assets)
%   (c) a panel with every observation missing has LogL=0 and LogLi=0
%   (d) the same panel given by Data.age and by Data.agej0 gives the same LogLi
%   (e) estimoptions.nindividualsbatch=7 gives the same LogLi as all individuals at once
%   (f) estimoptions.initialcondition='conditional' gives LogLi minus the first-age term LogLij(i,j0(i))
%   (g) runtime of a larger panel (informational)

Params.mewj=0.97.^(0:1:N_j-1); Params.mewj=Params.mewj/sum(Params.mewj);
jequaloneDist_z=jequaloneDist;
vfoptions_z=vfoptions;
simoptions_z=simoptions;
ReturnFn_z=ReturnFn;

for variant=1:2
clear FnsToEvaluate MeasurementError Data estimoptions
if variant==1
    vname='z';
    vfoptions=vfoptions_z; simoptions=simoptions_z; jequaloneDist=jequaloneDist_z; ReturnFn=ReturnFn_z;
    FnsToEvaluate.assets=@(aprime,a,z) a;
    FnsToEvaluate.cons=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) (1+r)*a+(agej<Jr)*w*kappa_j*z+(agej>=Jr)*pension-aprime;
    FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j,agej,Jr) (agej<Jr)*w*kappa_j*z;
    N_e=1;
else
    vname='ze';
    n_e=3;
    [e_grid,pi_e]=discretizeAR1_FarmerToda(0,0,0.1,n_e);
    pi_e=pi_e(1,:)';
    e_grid=exp(e_grid);
    vfoptions=vfoptions_z; simoptions=simoptions_z;
    vfoptions.n_e=n_e; vfoptions.e_grid=e_grid; vfoptions.pi_e=pi_e;
    simoptions.n_e=n_e; simoptions.e_grid=e_grid; simoptions.pi_e=pi_e;
    jequaloneDist=jequaloneDist_z.*reshape(gpuArray(pi_e),[1,1,n_e]);
    ReturnFn=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
    FnsToEvaluate.assets=@(aprime,a,z,e) a;
    FnsToEvaluate.cons=@(aprime,a,z,e,r,w,kappa_j,agej,Jr,pension) (1+r)*a+(agej<Jr)*w*kappa_j*z*e+(agej>=Jr)*pension-aprime;
    FnsToEvaluate.earnings=@(aprime,a,z,e,w,kappa_j,agej,Jr) (agej<Jr)*w*kappa_j*z*e;
    N_e=n_e;
end

[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);

N_a=prod(n_a);
N_z=prod(n_z);
N_s=N_a*N_z*N_e;
mewj=Params.(AgeWeightParamNames{1});
SD=gather(reshape(StationaryDist,[N_s,N_j]))./mewj; % each age of mass one
Vcons=gather(reshape(ValuesOnGrid.cons,[N_s,N_j]));
Vearn=gather(reshape(ValuesOnGrid.earnings,[N_s,N_j]));
Vassets=gather(reshape(ValuesOnGrid.assets,[N_s,N_j]));
agrid=gather(a_grid);

MeasurementError.cons.distribution='normal'; MeasurementError.cons.sigma='sig_c';
MeasurementError.earnings.distribution='t'; MeasurementError.earnings.sigma='sig_y'; MeasurementError.earnings.dof='nu';

%% (a) PopDist vs StationaryDist
Data.cons=NaN(3,1); Data.age=[1;2;3];
MEa.cons=MeasurementError.cons;
[~,~,~,PopDist]=LogLikelihood_FHorz(Data,MEa,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions);
temp=gather(reshape(PopDist,[N_s,N_j])).*mewj-gather(reshape(StationaryDist,[N_s,N_j]));
fprintf(['SolvedModel ',vname,' (a) PopDist times the age weights vs StationaryDist, this should be zero: %.3e \n'],CoreMaxAbs(temp))
fprintf(['SolvedModel ',vname,' (a) (discriminator) PopDist without the age weights vs StationaryDist is %.3e, and the mass of PopDist at each age minus one is at most %.3e (informational) \n'], ...
    max(abs(gather(reshape(PopDist,[N_s,N_j]))-gather(reshape(StationaryDist,[N_s,N_j]))),[],'all'),max(abs(sum(gather(reshape(PopDist,[N_s,N_j])),1)-1)))

%% (b) Cross-section vs the closed-form mixture
N=30;
clear Data
Data.age=mod((0:1:N-1)',N_j)+1;
Data.cons=zeros(N,1); Data.earnings=zeros(N,1); Data.assets=zeros(N,1);
for ii=1:N
    jj=Data.age(ii);
    mc=sum(SD(:,jj).*Vcons(:,jj)); sc=sqrt(sum(SD(:,jj).*(Vcons(:,jj)-mc).^2));
    Data.cons(ii)=mc+0.7*sc*sin(2*ii)+0.05*cos(ii);
    me=sum(SD(:,jj).*Vearn(:,jj));
    Data.earnings(ii)=me+0.1*sin(3*ii);
    asupport=find(sum(reshape(SD(:,jj),[N_a,N_z*N_e]),2)>0); % a points with mass at this age
    asupport=asupport(asupport<N_a);
    kk=asupport(mod(ii,numel(asupport))+1);
    Data.assets(ii)=agrid(kk)+(0.15+0.7*mod(ii,2))*(agrid(kk+1)-agrid(kk));
end
Data.cons(4)=NaN; Data.earnings(7)=NaN; Data.assets(9)=NaN; Data.cons(12)=NaN; Data.assets(12)=NaN;
MEb=MeasurementError;
MEb.assets.distribution='gridlottery';
[~,LogLi]=LogLikelihood_FHorz(Data,MEb,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions);
LogLi=gather(LogLi);
LogLi_cf=zeros(N,1);
sg=Params.sig_c; sy=Params.sig_y; nu=Params.nu;
for ii=1:N
    jj=Data.age(ii);
    f=ones(N_s,1);
    if ~isnan(Data.cons(ii))
        f=f.*exp(-0.5*((Data.cons(ii)-Vcons(:,jj))/sg).^2)/(sg*sqrt(2*pi));
    end
    if ~isnan(Data.earnings(ii))
        f=f.*gamma((nu+1)/2)/(gamma(nu/2)*sqrt(nu*pi)*sy).*(1+((Data.earnings(ii)-Vearn(:,jj))/sy).^2/nu).^(-(nu+1)/2);
    end
    if ~isnan(Data.assets(ii))
        kk=find(agrid<=Data.assets(ii),1,'last');
        wgt=(Data.assets(ii)-agrid(kk))/(agrid(kk+1)-agrid(kk));
        f=f.*((Vassets(:,jj)==agrid(kk))*(1-wgt)+(Vassets(:,jj)==agrid(kk+1))*wgt);
    end
    LogLi_cf(ii)=log(sum(SD(:,jj).*f));
end
fprintf(['SolvedModel ',vname,' (b) cross-section LogLi vs the closed-form mixture over StationaryDist, this should be zero: %.3e \n'],CoreMaxAbs(LogLi-LogLi_cf))
fprintf(['SolvedModel ',vname,' (b) LogLi ranges from %.4f to %.4f (informational) \n'],min(LogLi),max(LogLi))

%% A panel (T=4, consecutive ages)
N=40; T=4;
agej0=mod((0:1:N-1)',N_j-T+1)+1;
Age=agej0+(0:1:T-1);
clear Data
Data.cons=zeros(N,T); Data.earnings=zeros(N,T);
for ii=1:N
    for tt=1:T
        jj=Age(ii,tt);
        mc=sum(SD(:,jj).*Vcons(:,jj)); sc=sqrt(sum(SD(:,jj).*(Vcons(:,jj)-mc).^2));
        Data.cons(ii,tt)=mc+0.7*sc*sin(2*ii+5*tt);
        Data.earnings(ii,tt)=sum(SD(:,jj).*Vearn(:,jj))+0.1*sin(3*ii+tt);
    end
end
Data.cons(3,2)=NaN; Data.earnings(3,2)=NaN; Data.cons(5,1)=NaN; Data.earnings(8,4)=NaN; Data.cons(10,:)=NaN;
DataA=Data; DataA.age=Age;
DataB=Data; DataB.agej0=agej0;

%% (c) every observation missing
DataNaN=DataA; DataNaN.cons=NaN(N,T); DataNaN.earnings=NaN(N,T);
[LogL,LogLi]=LogLikelihood_FHorz(DataNaN,MeasurementError,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions);
fprintf(['SolvedModel ',vname,' (c) every observation missing: LogL, this should be zero: %.3e, LogLi, this should be zero: %.3e \n'],abs(gather(LogL)),CoreMaxAbs(gather(LogLi)))

%% (d) Data.age vs Data.agej0
[~,LogLiA,LogLijA]=LogLikelihood_FHorz(DataA,MeasurementError,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions);
[~,LogLiB]=LogLikelihood_FHorz(DataB,MeasurementError,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions);
LogLiA=gather(LogLiA); LogLijA=gather(LogLijA); LogLiB=gather(LogLiB);
fprintf(['SolvedModel ',vname,' (d) panel LogLi with Data.age vs with Data.agej0, this should be zero: %.3e \n'],CoreMaxAbs(LogLiA-LogLiB))
fprintf(['SolvedModel ',vname,' (d) panel LogL=%.4f; individual 10 (consumption never observed) has LogLi=%.4f, and the panel LogLi is not all equal: spread %.4f (informational) \n'],sum(LogLiA),LogLiA(10),max(LogLiA)-min(LogLiA))

%% (e) batches of individuals
estimoptions.nindividualsbatch=7;
[~,LogLiE]=LogLikelihood_FHorz(DataA,MeasurementError,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
fprintf(['SolvedModel ',vname,' (e) nindividualsbatch=7 vs all at once, this should be zero: %.3e \n'],CoreMaxAbs(gather(LogLiE)-LogLiA))

%% (f) conditional initial condition
clear estimoptions
estimoptions.initialcondition='conditional';
[~,LogLiF]=LogLikelihood_FHorz(DataA,MeasurementError,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
LogLiF=gather(LogLiF);
firstterm=LogLijA(sub2ind([N,N_j],(1:1:N)',agej0));
fprintf(['SolvedModel ',vname,' (f) conditional LogLi vs LogLi minus the first-age term, this should be zero: %.3e \n'],CoreMaxAbs(LogLiF-(LogLiA-firstterm)))
fprintf(['SolvedModel ',vname,' (f) (discriminator) the first-age terms are as large as %.4f in absolute value (informational) \n'],max(abs(firstterm)))

%% (g) runtime
N=2000; T=5;
agej0=mod((0:1:N-1)',N_j-T+1)+1;
clear Data
Data.agej0=agej0;
Data.cons=zeros(N,T); Data.earnings=zeros(N,T);
for tt=1:T
    jj=agej0+tt-1;
    for j=1:N_j
        ind=(jj==j);
        mc=sum(SD(:,j).*Vcons(:,j));
        Data.cons(ind,tt)=mc+0.3*sin((1:sum(ind))'+tt);
        Data.earnings(ind,tt)=sum(SD(:,j).*Vearn(:,j))+0.1*cos((1:sum(ind))'+tt);
    end
end
tic;
LogL=LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions);
LogL=gather(LogL);
temptime=toc;
fprintf(['SolvedModel ',vname,' (g) runtime: N=%i, T=%i, %i latent states: %.3f seconds, LogL=%.4f (informational) \n'],N,T,N_s,temptime,LogL)

output.(['LogLiA_',vname])=LogLiA;
end % variant

end
