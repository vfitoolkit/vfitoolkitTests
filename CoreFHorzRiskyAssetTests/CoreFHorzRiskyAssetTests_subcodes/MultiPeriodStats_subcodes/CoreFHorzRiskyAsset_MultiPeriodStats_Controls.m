function output=CoreFHorzRiskyAsset_MultiPeriodStats_Controls(n_d,d_grid,n_a,a_grid,n_z,N_j,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline)
% The option checks of the multi-period statistics for riskyasset, in EvalFnOnAgentDist_AutoCorrTransProbs_FHorz and the
% M>1 windows of EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz, EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz
% and LifeCycleProfiles_FHorz_Case1, checked as a set (every command):
%   positive control: riskyasset runs;
%   riskyasset together with experienceasset: refused ('only one of');
%   riskyasset without n_u, or without pi_u: refused (says n_u, u_grid and pi_u are needed);
% and in AutoCorrTransProbs (whose option checks are all at the top; in the other three the semiz check sits in the M>1 block,
% after the single-period statistics, which a call with a semiz but a Policy without one would already trip):
%   riskyasset with semiz but no refine_d: refused ('without setting simoptions.refine_d').
% Plus a PType cell: two identical types of the riskyasset model, EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType per type
% must equal the solo command (every numeric field; a NaN against a number below 1e-12 is the same answer, see
% CoreFHorzExpAsset_MultiPeriodStats_FakeExpAsset in the ExpAsset bank).
% (That residualasset and inheritanceasset are refused is checked in CoreFHorzExpAssetTests doPart(20); the numbers are checked
% in CoreFHorzRiskyAsset_MultiPeriodStats_BruteForce.)
% The model is the bank's riskyasset one with no d1, no a1, z and no e (riskyshare and savings, refine_d=[0,1,1]).

fprintf('\n========== CoreFHorzRiskyAsset_MultiPeriodStats_Controls ==========\n')

vfoptions=struct();
vfoptions.riskyasset=1;
vfoptions.refine_d=[0,1,1]; % no d1, 1 d2 (riskyshare), 1 d3 (savings)
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn; % (1+r)*(1-riskyshare)*savings+(1+r+u)*riskyshare*savings
vfoptions.n_u=vfoptionsbaseline.n_u;
vfoptions.u_grid=vfoptionsbaseline.u_grid;
vfoptions.pi_u=vfoptionsbaseline.pi_u;
simoptions=vfoptions;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;
simoptions.pi_z=pi_z;
ReturnFn=@(savings,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_noe_nosemiz(savings,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
F=struct();
F.assets=@(riskyshare,savings,a,z) a;
F.savings=@(riskyshare,savings,a,z) savings;
F.riskysavings=@(riskyshare,savings,a,z) riskyshare*savings;
jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(ceil(n_a/2),ceil(n_z/2))=1; % mid-grid wealth, so that the agents have savings to put in the risky asset from age 1

[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

cases={'riskyasset runs (positive control)','riskyasset and experienceasset both set','riskyasset without n_u','riskyasset without pi_u','riskyasset with semiz but no refine_d'};
expected={'','only one of','n_u, simoptions.u_grid and simoptions.pi_u','n_u, simoptions.u_grid and simoptions.pi_u','without setting simoptions.refine_d'};
cmds={'AutoCorrTransProbs_FHorz','CrossSectionCovarCorr_FHorz (M=2)','AgeConditionalStats_CrossSectionCovarCorr_FHorz (M=2)','LifeCycleProfiles_FHorz_Case1 (M=2)'};
nchecks=0;
for kk=1:length(cases)
    for cmd=1:length(cmds)
        if kk==5 && cmd>1 % the semiz check of the other three sits after their single-period statistics (see above)
            continue
        end
        so=simoptions;
        so.ntimeperiods=[1,2];
        if kk==2
            so.experienceasset=1;
        elseif kk==3
            so=rmfield(so,'n_u');
        elseif kk==4
            so=rmfield(so,'pi_u');
        elseif kk==5
            so.n_semiz=2;
            so=rmfield(so,'refine_d');
        end
        errmsg='';
        try
            if cmd==1
                EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,F,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,so);
            elseif cmd==2
                EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,F,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,so);
            elseif cmd==3
                EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,F,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,so);
            else
                LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,F,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,so);
            end
        catch ME
            errmsg=ME.message;
        end
        nchecks=nchecks+1;
        if kk==1
            bad=~isempty(errmsg);
            fprintf('MultiPeriodStats riskyasset controls: %s, %s, errored, this should be zero: %i \n',cmds{cmd},cases{kk},bad)
        else
            bad=~contains(errmsg,expected{kk});
            fprintf('MultiPeriodStats riskyasset controls: %s, %s, not refused with the expected error, this should be zero: %i \n',cmds{cmd},cases{kk},bad)
        end
        if bad && ~isempty(errmsg)
            fprintf('    the error was: %s \n',errmsg)
        end
    end
end
fprintf('MultiPeriodStats riskyasset controls: %i (command, case) pairs checked, of %i \n',nchecks,4*length(cmds)+1)

%% PType: two identical types of the riskyasset model, per type equal to the solo command
Names_i={'typeA','typeB'};
PTypeDistParamNames={'ptypemass'};
ParamsPT=Params;
ParamsPT.ptypemass=[0.5;0.5];
[~,PolicyPT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsPT,DiscountFactorParamNames,vfoptions);
StationaryDistPT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,PolicyPT,n_d,n_a,n_z,N_j,Names_i,pi_z,ParamsPT,simoptions);
soPT=simoptions; soPT.timehorizons=[1,2]; soPT.ntimeperiods=[1,2];
ACPT=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDistPT,PolicyPT,F,ParamsPT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,soPT);
ACsolo=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,F,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,soPT);
fprintf('MultiPeriodStats riskyasset PType: assets has zero StdDeviation at every age (the agents never move, the comparison would be vacuous), this should be zero: %i \n',all(gather(ACsolo.assets.StdDeviation)<1e-12 | isnan(gather(ACsolo.assets.StdDeviation))))
soloNames=fieldnames(ACsolo);
for ii=1:length(Names_i)
    maxdiff=0; worst='(none)'; nnanbad=0; nmismatch=0; nnumeric=0;
    for kk=1:length(soloNames)
        if ~isstruct(ACsolo.(soloNames{kk}))
            continue
        end
        if ~isfield(ACPT,soloNames{kk}) || ~isfield(ACPT.(soloNames{kk}),Names_i{ii})
            nmismatch=nmismatch+1;
            continue
        end
        X=ACPT.(soloNames{kk}).(Names_i{ii}); Y=ACsolo.(soloNames{kk});
        statNames=fieldnames(Y);
        for ss=1:length(statNames)
            if ~isfield(X,statNames{ss}) || ~(isnumeric(Y.(statNames{ss})) || islogical(Y.(statNames{ss})))
                nmismatch=nmismatch+~isfield(X,statNames{ss});
                continue
            end
            Xd=gather(double(X.(statNames{ss}))); Yd=gather(double(Y.(statNames{ss})));
            if ~isequal(size(Xd),size(Yd))
                nmismatch=nmismatch+1;
                continue
            end
            nnumeric=nnumeric+numel(Xd);
            nanbad=find(isnan(Xd(:))~=isnan(Yd(:))); % a NaN against double noise (below 1e-12) is the same answer
            for kk2=1:length(nanbad)
                defined=[Xd(nanbad(kk2)),Yd(nanbad(kk2))];
                defined=defined(~isnan(defined));
                nnanbad=nnanbad+(abs(defined)>1e-12);
                fprintf('    NaN in only one: %s.%s, entry %i: PType %.6e, solo %.6e \n',soloNames{kk},statNames{ss},nanbad(kk2),Xd(nanbad(kk2)),Yd(nanbad(kk2)))
            end
            both=~isnan(Xd(:)) & ~isnan(Yd(:));
            if any(both)
                dd=max(abs(Xd(both)-Yd(both)));
                if dd>maxdiff
                    maxdiff=dd;
                    worst=[soloNames{kk},'.',statNames{ss}];
                end
            end
        end
    end
    fprintf('MultiPeriodStats riskyasset PType (type %s): %i numbers compared, largest difference at %s \n',Names_i{ii},nnumeric,worst)
    fprintf('MultiPeriodStats riskyasset PType (type %s): AutoCorrTransProbs_FHorz_PType per type vs solo, this should be zero: %.3e \n',Names_i{ii},maxdiff)
    fprintf('MultiPeriodStats riskyasset PType (type %s): entries NaN in only one, this should be zero: %i \n',Names_i{ii},nnanbad)
    fprintf('MultiPeriodStats riskyasset PType (type %s): fields missing or of a different size, this should be zero: %i \n',Names_i{ii},nmismatch)
    fprintf('MultiPeriodStats riskyasset PType (type %s): nothing compared (a vacuous pass), this should be zero: %i \n',Names_i{ii},nnumeric==0)
end

output=struct();

end
