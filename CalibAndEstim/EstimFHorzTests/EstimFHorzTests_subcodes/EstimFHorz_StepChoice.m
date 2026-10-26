function output=EstimFHorz_StepChoice(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% EstimFHorzTests part 14: estimoptions.eedefault='auto' (the default) in EstimateLifeCycleModel_MethodOfMoments, which chooses the
% finite-difference step of each parameter: epsilon10000 (eedefault=3) where the up and down one-sided Jacobian columns agree there
% (estsummary.stepchoice.onesideddisagreement, ||J_up-J_down||_W/||J_up+J_down||_W, at most 0.1), and epsilon1000000 (eedefault=4)
% otherwise. The model moments are step functions of the parameters (the policy is on a grid), and a step inside their jumps gives a
% Jacobian that is mostly jumps.
% The targets are the moments of part 1 at the truth with the grid interpolation layer (ngridinterp=15); the standard errors at the
% truth (skipestimation=1) with eedefault=3, 4 and 'auto', on four model grids, from coarse to fine: no grid interpolation layer,
% ngridinterp=3, 15 (the targets') and 60. W the identity, Omega an arbitrary positive-definite diagonal, as part 1.
% Identities, on each grid:
%   (a) onesideddisagreement vs the test's own computation from doublechecks J_up, J_down and W; eechosen vs the rule applied to it;
%       eedefault=3 and 4 give eechosen all 3 and all 4
%   (b) J: each column is that of doublechecks at the parameter's chosen step; Sigma vs inv(J'WJ) J'W Omega W J inv(J'WJ) from it;
%       stepchoice.stderr vs sqrt(diag(Sigma)) of doublechecks
% Informational: on each grid, the step chosen, the one-sided disagreement at each step, the standard errors at 3, 4 and 'auto', and
% relative to those of 'auto' on the finest grid.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;

vfoptions=struct(); simoptions=struct();
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
simoptions_full=simoptions; simoptions_full.whichstats=ones(1,7);

%% Targets at the truth (as part 1)
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
TM=struct();
t=gather(LCPFull.assets.Mean); t(1)=NaN; t(t<1e-6)=NaN; TM.AgeConditionalStats.assets.Mean=t;
nomit=sum(isnan(t));
TM.AllStats.assets.Gini=gather(AllFull.assets.Gini);
TM.AllStats.earnings.Mean=gather(AllFull.earnings.Mean);
nT=(N_j-nomit)+2;
W=eye(nT);
Omega=diag(0.01*(1+(1:nT)'/nT));

EstimParamNames={'beta','pension'};
K=2;
estimoptions=struct();
estimoptions.verbose=0;
estimoptions.skipestimation=1;
epsnames={'epsilon1','epsilon100','epsilon10000','epsilon1000000'};
gridnames={'no grid interpolation layer','ngridinterp=3','ngridinterp=15','ngridinterp=60'};
eeruns={3,4,'auto'};
eerunnames={'3','4','auto'};

sd=zeros(K,3,4); % parameter, run (3, 4, auto), grid
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
        [~,~,es{ss}]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,esto,vfo,sio);
    end
    ea=es{3}; % the 'auto' run

    %% (a) the one-sided disagreement and the rule
    rhand=zeros(K,4);
    for ee=1:4
        dc=ea.doublechecks.(epsnames{ee});
        for pp=1:K
            dJ=dc.J_up(:,pp)-dc.J_down(:,pp); sJ=dc.J_up(:,pp)+dc.J_down(:,pp);
            rhand(pp,ee)=sqrt((dJ'*W*dJ)/(sJ'*W*sJ));
        end
    end
    temp=ea.stepchoice.onesideddisagreement-rhand;
    temp(isnan(ea.stepchoice.onesideddisagreement) & isnan(rhand))=0;
    eehand=3+(~(rhand(:,3)<=0.1));
    ec=ea.stepchoice.eechosen;
    fprintf('StepChoice (a) %s: onesideddisagreement vs ||J_up-J_down||_W/||J_up+J_down||_W from doublechecks, this should be zero: %.3e; eechosen vs the rule, this should be zero: %i; the fixed-step runs (eedefault three and four) keep their step, this should be zero: %i \n', ...
        gridnames{gg},CoreMaxAbs(temp),sum(ec(:)~=eehand),sum(es{1}.stepchoice.eechosen(:)~=3)+sum(es{2}.stepchoice.eechosen(:)~=4))

    %% (b) J, Sigma, stderr
    temp=0;
    for pp=1:K
        temp=temp+CoreMaxAbs(ea.variousmatrices.J(:,pp)-ea.doublechecks.(epsnames{ec(pp)}).J(:,pp));
    end
    J=ea.variousmatrices.J;
    SigmaHand=((J'*W*J)\(J'*W))*Omega*((J'*W*J)\(J'*W))';
    temp2=0;
    for ee=1:4
        temp3=sqrt(diag(ea.doublechecks.(epsnames{ee}).Sigma));
        temp4=ea.stepchoice.stderr(:,ee)-temp3;
        temp4(ea.stepchoice.stderr(:,ee)==temp3 | (isnan(ea.stepchoice.stderr(:,ee)) & isnan(temp3)))=0; % (a zero J column gives NaN or Inf on both sides)
        temp2=temp2+CoreMaxAbs(temp4);
    end
    fprintf('StepChoice (b) %s: J columns vs doublechecks at the chosen steps, this should be zero: %.3e; Sigma vs the sandwich from J, relative, this should be zero: %.3e; stderr vs sqrt(diag(Sigma)) of doublechecks, this should be zero: %.3e \n', ...
        gridnames{gg},temp,CoreMaxAbs(ea.variousmatrices.Sigma(:)-SigmaHand(:))/CoreMaxAbs(SigmaHand(:)),temp2)

    %% Informational
    fprintf('StepChoice %s: eechosen [beta, pension] = [%i, %i] (informational) \n',gridnames{gg},ec(1),ec(2))
    for pp=1:K
        fprintf('StepChoice %s, %s: one-sided disagreement at steps 1-4: %.3e %.3e %.3e %.3e (informational) \n',gridnames{gg},EstimParamNames{pp},rhand(pp,:))
    end
    for ss=1:3
        sd(:,ss,gg)=sqrt(diag(es{ss}.variousmatrices.Sigma));
        fprintf('StepChoice %s, eedefault=%s: s.e. beta %.3e pension %.3e (informational) \n',gridnames{gg},eerunnames{ss},sd(1,ss,gg),sd(2,ss,gg))
    end
end

%% Informational: standard errors relative to those of 'auto' on the finest grid
for gg=1:4
    for ss=1:3
        fprintf('StepChoice %s, eedefault=%s: s.e. relative to auto on ngridinterp=60, beta %.3f pension %.3f (informational) \n',gridnames{gg},eerunnames{ss},sd(1,ss,gg)/sd(1,3,4),sd(2,ss,gg)/sd(2,3,4))
    end
end

output.sd=sd;
end
