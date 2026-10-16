function output=EstimFHorzMLE_StepChoice(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions)
% doPart(10): estimoptions.eedefault='auto' (the default), which chooses the finite-difference step of each parameter: epsilon10000
% (eedefault=3) where the up and down one-sided scores agree there (estsummary.stepchoice.onesideddisagreement at most 0.1), and
% epsilon1000000 (eedefault=4) otherwise. The likelihood is a step function of the parameters (the policy is on a grid), and a step
% inside its jumps inflates the OPG and gives too-small standard errors.
% One panel (N=1000) simulated by hand at the truth with the bank's grid (ngridinterp=15), and the standard errors at the truth
% (skipestimation=1) with eedefault=3, 4 and 'auto', on four estimation grids, from coarse to fine: no grid interpolation layer,
% ngridinterp=3, 15 (the DGP's) and 60.
% Identities, on each grid:
%   (a) onesideddisagreement vs the test's own one-sided scores from doublechecks; eechosen vs the rule applied to them; eedefault=3
%       and 4 give eechosen all 3 and all 4; the doublechecks are the same whatever eedefault
%   (b) Scores: each column is that of doublechecks at the parameter's chosen step; stderr_OPG vs sqrt(diag(Sigma_OPG)) of doublechecks
%   (c) the Hessian: each diagonal entry equals that of the eedefault run at the parameter's chosen step; the off-diagonal entry equals
%       that run's when both parameters chose the same step, and otherwise its formula from the reported points and steps
% Informational (how 'auto' performs; there is no Monte Carlo, so no coverage): on each grid, the step chosen, the one-sided
% disagreement and the share of individuals with a score of exactly zero at each step, the OPG and Hessian standard errors at 3, 4
% and 'auto', and the information-matrix ratio diag(OPG)./diag(-Hessian) (near one when both are reliable, as the model is correctly
% specified); and the standard errors relative to those of 'auto' on the finest grid (where the jumps are smallest).

Params.sig_c=0.1; Params.sig_y=0.05; Params.sig_a=0.2; % the truth (measurement error)
FnsToEvaluate.cons=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) (1+r)*a+(agej<Jr)*w*kappa_j*z+(agej>=Jr)*pension-aprime;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j,agej,Jr) (agej<Jr)*w*kappa_j*z;
FnsToEvaluate.assets=@(aprime,a,z) a;
MeasurementError.cons.distribution='normal'; MeasurementError.cons.sigma='sig_c';
MeasurementError.earnings.distribution='normal'; MeasurementError.earnings.sigma='sig_y';
MeasurementError.assets.distribution='normal'; MeasurementError.assets.sigma='sig_a';

[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
N=1000;
Data=EstimFHorzMLE_SimulatePanel(Policy,n_a,n_z,N_j,a_grid,z_grid,pi_z,jequaloneDist,Params,vfoptions.ngridinterp,N,20261016);

EstimParamNames={'beta','sig_c'};
K=2;
estimoptions.verbose=0;
estimoptions.skipestimation=1;
estimoptions.constrainpositive={'sig_c'};
epsnames={'epsilon1','epsilon100','epsilon10000','epsilon1000000'};
gridnames={'no grid interpolation layer','ngridinterp=3','ngridinterp=15','ngridinterp=60'};
eeruns={3,4,'auto'};
eerunnames={'3','4','auto'};

sd_OPG=zeros(K,3,4); % parameter, run (3, 4, auto), grid
sd_H=zeros(K,3,4);
for gg=1:4
    vfo=vfoptions; sio=simoptions;
    if gg==1
        vfo.gridinterplayer=0; sio.gridinterplayer=0;
    else
        ngi=[0,3,15,60];
        vfo.ngridinterp=ngi(gg); sio.ngridinterp=ngi(gg);
    end
    es=cell(3,1);
    for ss=1:3
        esto=estimoptions; esto.eedefault=eeruns{ss};
        [~,~,es{ss}]=EstimateLifeCycleModel_MaximumLikelihood(EstimParamNames,Data,MeasurementError,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,[],FnsToEvaluate,esto,vfo,sio);
    end
    ea=es{3}; % the 'auto' run
    LogLi=ea.doublechecks.LogLi;
    x=ea.doublechecks.epsilon10000.estimparamsvec(:);

    %% (a) the one-sided disagreement and the rule
    rhand=zeros(K,4); zeroshare=zeros(K,4);
    for ee=1:4
        dc=ea.doublechecks.(epsnames{ee});
        Su=(dc.LogLi_up-LogLi)./((dc.estimparamsvecup(:)-x)');
        Sd=(LogLi-dc.LogLi_down)./((x-dc.estimparamsvecdown(:))');
        rhand(:,ee)=sqrt(sum((Su-Sd).^2,1)./sum((Su+Sd).^2,1))';
        zeroshare(:,ee)=mean(dc.LogLi_up==dc.LogLi_down,1)';
    end
    temp=ea.stepchoice.onesideddisagreement-rhand;
    temp(isnan(ea.stepchoice.onesideddisagreement) & isnan(rhand))=0;
    eehand=3+(~(rhand(:,3)<=0.1));
    fprintf('StepChoice (a) %s: onesideddisagreement vs the one-sided scores from doublechecks, this should be zero: %.3e; eechosen vs the rule, this should be zero: %i \n',gridnames{gg},CoreMaxAbs(temp),sum(ea.stepchoice.eechosen(:)~=eehand))
    temp=0;
    for ee=1:4
        temp=temp+CoreMaxAbs(es{1}.doublechecks.(epsnames{ee}).scores-ea.doublechecks.(epsnames{ee}).scores)+CoreMaxAbs(es{2}.doublechecks.(epsnames{ee}).scores-ea.doublechecks.(epsnames{ee}).scores);
    end
    fprintf('StepChoice (a) %s: the fixed-step runs (eedefault three and four) keep their step, this should be zero: %i; doublechecks scores the same whatever eedefault, this should be zero: %.3e \n',gridnames{gg},sum(es{1}.stepchoice.eechosen(:)~=3)+sum(es{2}.stepchoice.eechosen(:)~=4),temp)

    %% (b) the scores and stderr_OPG
    temp=0; temp2=0;
    for pp=1:K
        temp=temp+CoreMaxAbs(ea.variousmatrices.Scores(:,pp)-ea.doublechecks.(epsnames{ea.stepchoice.eechosen(pp)}).scores(:,pp));
    end
    for ee=1:4
        temp3=sqrt(diag(ea.doublechecks.(epsnames{ee}).Sigma_OPG));
        temp4=ea.stepchoice.stderr_OPG(:,ee)-temp3;
        temp4(ea.stepchoice.stderr_OPG(:,ee)==temp3 | (isnan(ea.stepchoice.stderr_OPG(:,ee)) & isnan(temp3)))=0; % (all-zero scores give NaN or Inf on both sides)
        temp2=temp2+CoreMaxAbs(temp4);
    end
    fprintf('StepChoice (b) %s: Scores columns vs doublechecks at the chosen steps, this should be zero: %.3e; stderr_OPG vs sqrt(diag(Sigma_OPG)) of doublechecks, this should be zero: %.3e \n',gridnames{gg},temp,temp2)

    %% (c) the Hessian
    H=ea.variousmatrices.Hessian;
    ec=ea.stepchoice.eechosen;
    temp=0;
    for pp=1:K
        temp=temp+abs(H(pp,pp)-es{ec(pp)-2}.variousmatrices.Hessian(pp,pp));
    end
    hs=ea.doublechecks.hessian;
    if ec(1)==ec(2)
        temp2=abs(H(1,2)-es{ec(1)-2}.variousmatrices.Hessian(1,2));
    else
        temp2=abs(H(1,2)-(hs.LogL_upup(1,2)-hs.LogL_up(1)-hs.LogL_up(2)+2*hs.LogL-hs.LogL_down(1)-hs.LogL_down(2)+hs.LogL_downdown(1,2))/(2*hs.hup(1)*hs.hup(2)));
    end
    temp3=0;
    for pp=1:K
        temp3=temp3+abs(hs.hup(pp)-(ea.doublechecks.(epsnames{ec(pp)}).estimparamsvecup(pp)-x(pp)));
    end
    fprintf('StepChoice (c) %s: Hessian diagonal vs the eedefault run at each chosen step, this should be zero: %.3e; off-diagonal, this should be zero: %.3e; the Hessian steps are the chosen ones, this should be zero: %.3e \n',gridnames{gg},temp,temp2,temp3)

    %% Informational: how the steps perform
    fprintf('StepChoice %s: eechosen [beta, sig_c] = [%i, %i] (informational) \n',gridnames{gg},ec(1),ec(2))
    for pp=1:K
        fprintf('StepChoice %s, %s: one-sided disagreement at steps 1-4: %.3e %.3e %.3e %.3e; share of zero scores: %.3f %.3f %.3f %.3f (informational) \n',gridnames{gg},EstimParamNames{pp},rhand(pp,:),zeroshare(pp,:))
    end
    for ss=1:3
        dOPG=diag(es{ss}.variousmatrices.Sigma_OPG);
        dH=diag(es{ss}.variousmatrices.Sigma_Hessian);
        sdH=sqrt(abs(dH)); sdH(dH<0)=NaN; % a Hessian that is not negative definite has no standard errors
        sd_OPG(:,ss,gg)=sqrt(dOPG);
        sd_H(:,ss,gg)=sdH;
        IMratio=diag(es{ss}.variousmatrices.OPG)./diag(-es{ss}.variousmatrices.Hessian);
        fprintf('StepChoice %s, eedefault=%s: s.e. (OPG) beta %.3e sig_c %.3e; s.e. (Hessian) beta %.3e sig_c %.3e; diag(OPG)./diag(-Hessian), near one when both are reliable: beta %.3f sig_c %.3f (informational) \n',gridnames{gg},eerunnames{ss},sd_OPG(1,ss,gg),sd_OPG(2,ss,gg),sd_H(1,ss,gg),sd_H(2,ss,gg),IMratio(1),IMratio(2))
    end
end

%% Informational: standard errors relative to those of 'auto' on the finest grid
for gg=1:4
    for ss=1:3
        fprintf('StepChoice %s, eedefault=%s: s.e. relative to auto on ngridinterp=60, (OPG) beta %.3f sig_c %.3f; (Hessian) beta %.3f sig_c %.3f (informational) \n',gridnames{gg},eerunnames{ss},sd_OPG(1,ss,gg)/sd_OPG(1,3,4),sd_OPG(2,ss,gg)/sd_OPG(2,3,4),sd_H(1,ss,gg)/sd_OPG(1,3,4),sd_H(2,ss,gg)/sd_OPG(2,3,4))
    end
end

output.sd_OPG=sd_OPG;
output.sd_H=sd_H;
end
