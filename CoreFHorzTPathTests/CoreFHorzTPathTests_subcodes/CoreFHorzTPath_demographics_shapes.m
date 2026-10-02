function output=CoreFHorzTPath_demographics_shapes(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,n_d2_semiz,d2_grid_semiz,Params,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline)
% (10) of the demographics checks: the age weights and the newborn distribution are reshaped in a
% separate branch for each combination of shocks (and each of fastOLG=0/1), which is where these go
% wrong. So checks (1), (2) and (6) of CoreFHorzTPath_demographics are repeated on a model without z,
% one with z and e, and one with a semi-exogenous state. See that file for what each check is.
% (6) is not done with the semi-exogenous state: jequalOneDist as a T-path is not implemented there.

DiscountFactorParamNames={'beta'};

%% The same demographic transition as CoreFHorzTPath_demographics
sj_initial=[linspace(0.99,0.9,N_j-1),0]';
sj_final=[linspace(0.995,0.95,N_j-1),0]';
n_initial=0.01;
n_final=0;
ramp=min((0:T-1)/(T/2),1);
sjPath=sj_initial+(sj_final-sj_initial).*ramp;
nPath=n_initial+(n_final-n_initial)*ramp;
massPath=zeros(N_j,T);
massPath(1,1)=1;
for jj=2:N_j
    massPath(jj,1)=massPath(jj-1,1)*sj_initial(jj-1)/(1+n_initial);
end
for tt=2:T
    massPath(1,tt)=massPath(1,tt-1)*(1+nPath(tt));
    massPath(2:N_j,tt)=massPath(1:N_j-1,tt-1).*sjPath(1:N_j-1,tt-1);
end
mewjPath=massPath./sum(massPath,1);
Params.mewj=mewjPath(:,1)';

ParamPath1=ParamPath;
ParamPath1.mewj=mewjPath;
ParamPath1O=ParamPath;
ParamPath1O.mewjratio=mewjPath./mewjPath(:,1);
ParamPath2=ParamPath;
ParamPath2.mewj=repmat(mewjPath(:,1),1,T);
Params1O=Params;
Params1O.mewjratio=ones(1,N_j);

ia_T=1+floor(20*(0:T-1)/T); % newborns enter at grid points 1 to 20 over the path

transpathoptionsE=transpathoptionsbaseline;
transpathoptionsE.maxiter=1;
transpathoptionsE.GEnewprice=3;
transpathoptionsE.GEnewprice3.howtoupdate={'g','r',0,0.1};
transpathoptionsE.verbose=0;
GeneralEqmEqnsAE.g=@(assets,earnings) assets+10*earnings;
GeneralEqmEqns6.g=@(assets1) assets1;
fastOLGopts=[1,0];


%% Model without z
vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];
n_z0=0; z_grid0=[]; pi_z0=[];
jequaloneDist=zeros(n_a,1,'gpuArray');
jequaloneDist(1)=1;
ReturnFn=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate=struct();
FnsToEvaluate.assets=@(aprime,a) a;
FnsToEvaluate.earnings=@(aprime,a,w,kappa_j) w*kappa_j;
FnsToEvaluate1O=struct();
FnsToEvaluate1O.assets=@(aprime,a,mewjratio) a*mewjratio;
FnsToEvaluate1O.earnings=@(aprime,a,w,kappa_j,mewjratio) w*kappa_j*mewjratio;
[V_final,Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z0,N_j,d_grid,a_grid,z_grid0,pi_z0,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_final,n_d,n_a,n_z0,N_j,pi_z0,Params,simoptions);
for ff=1:length(fastOLGopts)
    transpathoptionsE.fastOLG=fastOLGopts(ff);
    [~,GEcondnPath1]=TransitionPath_Case1_FHorz(PricePath, ParamPath1, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z0, N_j, d_grid,a_grid,z_grid0, pi_z0, ReturnFn, FnsToEvaluate, GeneralEqmEqnsAE, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    [~,GEcondnPath1O]=TransitionPath_Case1_FHorz(PricePath, ParamPath1O, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z0, N_j, d_grid,a_grid,z_grid0, pi_z0, ReturnFn, FnsToEvaluate1O, GeneralEqmEqnsAE, Params1O, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('Demographics, no z (fastOLG=%i): (1) mewj on the ParamPath vs period-1 mewj with the ratio in the FnsToEvaluate, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(GEcondnPath1.g-GEcondnPath1O.g)))
    [~,GEcondnPath2]=TransitionPath_Case1_FHorz(PricePath, ParamPath2, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z0, N_j, d_grid,a_grid,z_grid0, pi_z0, ReturnFn, FnsToEvaluate, GeneralEqmEqnsAE, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    [~,GEcondnPath2O]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z0, N_j, d_grid,a_grid,z_grid0, pi_z0, ReturnFn, FnsToEvaluate, GeneralEqmEqnsAE, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('Demographics, no z (fastOLG=%i): (2) constant mewj on the ParamPath vs not on it, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(GEcondnPath2.g-GEcondnPath2O.g)))
    jequaloneDist_T=zeros([n_a,T],'gpuArray');
    for tt=1:T
        jequaloneDist_T(ia_T(tt),tt)=1;
    end
    FnsToEvaluate6=struct();
    FnsToEvaluate6.assets1=@(aprime,a,agej) a*(agej==1);
    [~,GEcondnPath6]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist_T, n_d, n_a, n_z0, N_j, d_grid,a_grid,z_grid0, pi_z0, ReturnFn, FnsToEvaluate6, GeneralEqmEqns6, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('Demographics, no z (fastOLG=%i): (6) jequaloneDist as a T-path, age-1 assets in periods 2 to T-1 vs by hand, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(gather(GEcondnPath6.g(2:T-1))-Params.mewj(1)*gather(a_grid(ia_T(2:T-1)))')))
end


%% Model with z and e
vfoptions=struct();
simoptions=struct();
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
simoptions.pi_e=simoptionsbaseline.pi_e;
n_e=vfoptions.n_e;
jequaloneDist=zeros([n_a,n_z,n_e],'gpuArray');
jequaloneDist(1,ceil(n_z/2),ceil(n_e/2))=1;
ReturnFn=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate=struct();
FnsToEvaluate.assets=@(aprime,a,z,e) a;
FnsToEvaluate.earnings=@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;
FnsToEvaluate1O=struct();
FnsToEvaluate1O.assets=@(aprime,a,z,e,mewjratio) a*mewjratio;
FnsToEvaluate1O.earnings=@(aprime,a,z,e,w,kappa_j,mewjratio) w*kappa_j*z*e*mewjratio;
[V_final,Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_final,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
for ff=1:length(fastOLGopts)
    transpathoptionsE.fastOLG=fastOLGopts(ff);
    [~,GEcondnPath1]=TransitionPath_Case1_FHorz(PricePath, ParamPath1, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsAE, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    [~,GEcondnPath1O]=TransitionPath_Case1_FHorz(PricePath, ParamPath1O, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate1O, GeneralEqmEqnsAE, Params1O, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('Demographics, z and e (fastOLG=%i): (1) mewj on the ParamPath vs period-1 mewj with the ratio in the FnsToEvaluate, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(GEcondnPath1.g-GEcondnPath1O.g)))
    [~,GEcondnPath2]=TransitionPath_Case1_FHorz(PricePath, ParamPath2, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsAE, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    [~,GEcondnPath2O]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsAE, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('Demographics, z and e (fastOLG=%i): (2) constant mewj on the ParamPath vs not on it, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(GEcondnPath2.g-GEcondnPath2O.g)))
    jequaloneDist_T=zeros([n_a,n_z,n_e,T],'gpuArray');
    for tt=1:T
        jequaloneDist_T(ia_T(tt),ceil(n_z/2),ceil(n_e/2),tt)=1;
    end
    FnsToEvaluate6=struct();
    FnsToEvaluate6.assets1=@(aprime,a,z,e,agej) a*(agej==1);
    [~,GEcondnPath6]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist_T, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate6, GeneralEqmEqns6, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('Demographics, z and e (fastOLG=%i): (6) jequaloneDist as a T-path, age-1 assets in periods 2 to T-1 vs by hand, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(gather(GEcondnPath6.g(2:T-1))-Params.mewj(1)*gather(a_grid(ia_T(2:T-1)))')))
end


%% Model with a semi-exogenous state (binary d2 drives it)
vfoptions=struct();
simoptions=struct();
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptions.d_grid=d2_grid_semiz;
n_d=n_d2_semiz;
d_grid=d2_grid_semiz;
jequaloneDist=zeros([n_a,vfoptions.n_semiz,n_z],'gpuArray');
jequaloneDist(1,ceil(vfoptions.n_semiz/2),ceil(n_z/2))=1;
ReturnFn=@(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
FnsToEvaluate=struct();
FnsToEvaluate.assets=@(d2,aprime,a,semiz,z) a;
FnsToEvaluate.earnings=@(d2,aprime,a,semiz,z,w,kappa_j) w*kappa_j*z*semiz;
FnsToEvaluate1O=struct();
FnsToEvaluate1O.assets=@(d2,aprime,a,semiz,z,mewjratio) a*mewjratio;
FnsToEvaluate1O.earnings=@(d2,aprime,a,semiz,z,w,kappa_j,mewjratio) w*kappa_j*z*semiz*mewjratio;
[V_final,Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_final,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
for ff=1:length(fastOLGopts)
    transpathoptionsE.fastOLG=fastOLGopts(ff);
    [~,GEcondnPath1]=TransitionPath_Case1_FHorz(PricePath, ParamPath1, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsAE, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    [~,GEcondnPath1O]=TransitionPath_Case1_FHorz(PricePath, ParamPath1O, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate1O, GeneralEqmEqnsAE, Params1O, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('Demographics, semiz (fastOLG=%i): (1) mewj on the ParamPath vs period-1 mewj with the ratio in the FnsToEvaluate, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(GEcondnPath1.g-GEcondnPath1O.g)))
    [~,GEcondnPath2]=TransitionPath_Case1_FHorz(PricePath, ParamPath2, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsAE, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    [~,GEcondnPath2O]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsAE, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('Demographics, semiz (fastOLG=%i): (2) constant mewj on the ParamPath vs not on it, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(GEcondnPath2.g-GEcondnPath2O.g)))
end

output=struct();

end
