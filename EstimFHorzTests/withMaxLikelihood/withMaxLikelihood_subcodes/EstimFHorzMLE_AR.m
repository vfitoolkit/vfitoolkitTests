function output=EstimFHorzMLE_AR(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions)
% doPart(7): AR(1) measurement error on the solved EstimFHorzTests model.
%   (a) PopDist with an AR(1) component equals PopDist without (the AR state is summed out, and must not change the population)
%   (b) a cross-section equals the closed-form mixture over StationaryDist and the stationary distribution of the AR(1) component
%   (c) a panel with every observation missing has LogL=0
%   (d) EstimateLifeCycleModel_MaximumLikelihood with rho_c and sigAR_c estimated, at the truth (skipestimation=1): LogL and the scores
%       of rho_c and sigAR_c vs direct calls at the reported up/down parameter vectors
%   (d2) the warning when the iid s.d. is smaller than the AR(1) grid spacing, and none when it is larger
%   (e) refusals

Params.sig_c=0.1; Params.sig_y=0.05; Params.sig_a=0.2;
Params.rho_c=0.7; Params.sigAR_c=0.05;
nAR=5;
FnsToEvaluate.cons=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) (1+r)*a+(agej<Jr)*w*kappa_j*z+(agej>=Jr)*pension-aprime;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j,agej,Jr) (agej<Jr)*w*kappa_j*z;
FnsToEvaluate.assets=@(aprime,a,z) a;
MeasurementError.cons.distribution='normal'; MeasurementError.cons.sigma='sig_c'; MeasurementError.cons.rho='rho_c'; MeasurementError.cons.sigmaAR='sigAR_c';
MeasurementError.earnings.distribution='normal'; MeasurementError.earnings.sigma='sig_y';
MeasurementError.assets.distribution='normal'; MeasurementError.assets.sigma='sig_a';
estimoptions.nARpoints=nAR;

[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
Params_sd=Params; Params_sd.mewj=ones(1,N_j)/N_j;
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,{'mewj'},Policy,n_d,n_a,n_z,N_j,pi_z,Params_sd,simoptions);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
N_s=prod(n_a)*prod(n_z);
SD=gather(reshape(StationaryDist,[N_s,N_j]))*N_j; % each age of mass one
Vcons=gather(reshape(ValuesOnGrid.cons,[N_s,N_j]));

%% (a) PopDist with and without the AR(1) component
Data.cons=NaN(3,1); Data.age=[1;2;3];
ME_noAR.cons=rmfield(MeasurementError.cons,{'rho','sigmaAR'});
MEa.cons=MeasurementError.cons;
[~,~,~,PopAR]=LogLikelihood_FHorz(Data,MEa,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
[~,~,~,PopNoAR]=LogLikelihood_FHorz(Data,ME_noAR,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
fprintf('AR (a) PopDist with an AR(1) component vs without, this should be zero: %.3e; vs StationaryDist, this should be zero: %.3e \n',CoreMaxAbs(gather(PopAR(:))-gather(PopNoAR(:))),CoreMaxAbs(reshape(gather(PopAR),[N_s,N_j])-SD))

%% (b) cross-section vs closed form
N=25;
clear Data
Data.age=mod((0:1:N-1)',N_j)+1;
Data.cons=zeros(N,1);
for ii=1:N
    jj=Data.age(ii);
    mc=sum(SD(:,jj).*Vcons(:,jj)); sc=sqrt(sum(SD(:,jj).*(Vcons(:,jj)-mc).^2));
    Data.cons(ii)=mc+0.7*sc*sin(2*ii);
end
[~,LogLi]=LogLikelihood_FHorz(Data,MEa,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
[ug,pu]=discretizeAR1_Rouwenhorst(0,Params.rho_c,Params.sigAR_c,nAR);
ug=gather(ug); pu=gather(pu);
[V,Dg]=eig(pu'); [~,kmax]=max(real(diag(Dg))); su=real(V(:,kmax)); su=su/sum(su);
LogLi_cf=zeros(N,1);
for ii=1:N
    jj=Data.age(ii);
    f=zeros(N_s,1);
    for uu=1:nAR
        f=f+su(uu)*exp(-0.5*((Data.cons(ii)-Vcons(:,jj)-ug(uu))/Params.sig_c).^2)/(Params.sig_c*sqrt(2*pi));
    end
    LogLi_cf(ii)=log(sum(SD(:,jj).*f));
end
fprintf('AR (b) cross-section LogLi vs the closed-form mixture over StationaryDist and the stationary AR(1) distribution, this should be zero: %.3e \n',CoreMaxAbs(gather(LogLi)-LogLi_cf))

%% (c) every observation missing
clear Data
Data.agej0=ones(10,1); Data.cons=NaN(10,5); Data.earnings=NaN(10,5); Data.assets=NaN(10,5);
[LogL_c,LogLi_c]=LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
fprintf('AR (c) every observation missing: LogL, this should be zero: %.3e; LogLi, this should be zero: %.3e \n',abs(gather(LogL_c)),CoreMaxAbs(gather(LogLi_c)))

%% (d) the estimator at the truth with rho_c and sigAR_c estimated
Data=EstimFHorzMLE_SimulatePanel(Policy,n_a,n_z,N_j,a_grid,z_grid,pi_z,jequaloneDist,Params,vfoptions.ngridinterp,200,20261020);
EstimParamNames={'rho_c','sigAR_c'};
estimoptions_d=estimoptions;
estimoptions_d.verbose=0;
estimoptions_d.skipestimation=1;
[~,~,es]=EstimateLifeCycleModel_MaximumLikelihood(EstimParamNames,Data,MeasurementError,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,[],FnsToEvaluate,estimoptions_d,vfoptions,simoptions);
dc=es.doublechecks.epsilon10000;
[LogLhand,LogLihand]=LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
LogLi_up=zeros(size(Data.cons,1),2); LogLi_down=LogLi_up;
for pp=1:2
    P=Params; P.rho_c=dc.estimparamsvec(1); P.sigAR_c=dc.estimparamsvec(2);
    P.(EstimParamNames{pp})=dc.estimparamsvecup(pp);
    [~,temp]=LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,Policy,FnsToEvaluate,P,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
    LogLi_up(:,pp)=gather(temp);
    P.(EstimParamNames{pp})=dc.estimparamsvecdown(pp);
    [~,temp]=LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,Policy,FnsToEvaluate,P,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
    LogLi_down(:,pp)=gather(temp);
end
Shand=(LogLi_up-LogLi_down)./((dc.estimparamsvecup-dc.estimparamsvecdown)');
fprintf('AR (d) estimator at the truth: LogL vs direct, this should be zero: %.3e; scores of rho_c and sigAR_c vs centered differences of direct calls, this should be zero: %.3e (relative to the largest score %.3e) \n',abs(es.LogL-gather(LogLhand)),CoreMaxAbs(es.variousmatrices.Scores-Shand),max(abs(Shand(:))))
fprintf('AR (d) standard errors (OPG): rho_c %.4f, sigAR_c %.4f (informational) \n',es.EstimParamsStdDev.rho_c,es.EstimParamsStdDev.sigAR_c)
fprintf('AR (d) reported AR grid spacing vs the spacing of the Rouwenhorst grid, this should be zero: %.3e; sigmaovergridspacing vs sig_c over it, this should be zero: %.3e (ratio %.3f, informational) \n',abs(es.measurementerror.cons.ARgridspacing-(ug(2)-ug(1))),abs(es.measurementerror.cons.sigmaovergridspacing-Params.sig_c/(ug(2)-ug(1))),es.measurementerror.cons.sigmaovergridspacing)

%% (d2) the warning when the iid s.d. is smaller than the AR grid spacing
lastwarn('');
P=Params; P.sig_c=0.01;
LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,Policy,FnsToEvaluate,P,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
fprintf('AR (d2) iid s.d. below the AR grid spacing warns, this should be one: %i \n',double(contains(lastwarn,'smaller than the spacing of the AR(1) grid')))
lastwarn('');
LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
fprintf('AR (d2) iid s.d. above the AR grid spacing does not warn, this should be zero: %i \n',double(contains(lastwarn,'smaller than the spacing of the AR(1) grid')))

%% (e) refusals
errcases={'AR component with exact','rho not a scalar','rho of one','missing sigmaAR','even nARpoints','nARpoints not set'};
errmsgs={'only possible with the normal and t distributions','must be a scalar','strictly between','needs both .rho and .sigmaAR','must be an odd number','you must set estimoptions.nARpoints'};
clear Data
Data.agej0=ones(3,1); Data.cons=ones(3,2);
for ee=1:length(errcases)
    ME=MEa; P=Params; esto=estimoptions;
    switch ee
        case 1
            ME.cons.distribution='exact';
        case 2
            P.rho_c=0.5*ones(1,N_j);
        case 3
            P.rho_c=1;
        case 4
            ME.cons=rmfield(ME.cons,'sigmaAR');
        case 5
            esto.nARpoints=4;
        case 6
            esto=rmfield(esto,'nARpoints');
    end
    errok=0;
    try
        LogLikelihood_FHorz(Data,ME,jequaloneDist,Policy,FnsToEvaluate,P,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions,esto);
    catch err
        errok=contains(err.message,errmsgs{ee});
        if ~errok
            fprintf('   (the error was: %s) \n',err.message);
        end
    end
    fprintf('AR (e) %s errors with the intended message, this should be zero: %i \n',errcases{ee},1-errok);
end

output.es=es;
end
