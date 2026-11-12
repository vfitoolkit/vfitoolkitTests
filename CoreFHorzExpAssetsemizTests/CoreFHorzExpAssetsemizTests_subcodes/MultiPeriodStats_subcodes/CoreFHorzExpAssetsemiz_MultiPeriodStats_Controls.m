function output=CoreFHorzExpAssetsemiz_MultiPeriodStats_Controls(n_d2,d2_grid,n_d3,d3_grid,n_a,a_grid,n_z,N_j,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline)
% The option checks of the multi-period statistics for experienceassetsemiz, in EvalFnOnAgentDist_AutoCorrTransProbs_FHorz and the
% M>1 windows of EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz, EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz
% and LifeCycleProfiles_FHorz_Case1, checked as a set (every command):
%   positive control: experienceassetsemiz runs;
%   experienceassetsemiz together with experienceasset, or with experienceassetu: refused ('only one of');
% and in AutoCorrTransProbs (whose option checks are all at the top; in the other three the shock check sits in the M>1 block, after
% the single-period statistics, which a call without the shock and a Policy over it would already trip):
%   experienceassetsemiz without the shock it depends on: refused (the message, 'experienceassete needs ... e, experienceassetze needs both z and e, ...', names what each needs).
% Plus a PType cell: two identical types of the experienceassetsemiz model, EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType per type
% must equal the solo command (every numeric field; a NaN against a number below 1e-12 is the same answer, see
% CoreFHorzExpAsset_MultiPeriodStats_FakeExpAsset in the ExpAsset bank).
% (That the other unhandled assets are refused is checked in CoreFHorzExpAssetTests doPart(20); the numbers are checked in
% CoreFHorzExpAssetsemiz_MultiPeriodStats_BruteForce.)
% The model is the bank's experienceassetsemiz one with no d1, a1, z and e and semiz (and its decision d3).

fprintf('\n========== CoreFHorzExpAssetsemiz_MultiPeriodStats_Controls ==========\n')

n_d=[n_d2,n_d3]; % d3 is the semi-exogenous decision
d_grid=[d2_grid;d3_grid];
vfoptions=struct();
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
vfoptions.experienceassetsemiz=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn; % aprimeFn of the bank
simoptions=vfoptions;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;
simoptions.pi_z=pi_z;
ReturnFn=@(d2,d3,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetsemiz_nod1_z_e(d2,d3,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
F.assets=@(d2,d3,a1prime,a1,a2,semiz,z,e) a1;
F.human=@(d2,d3,a1prime,a1,a2,semiz,z,e) a2;
F.earnings=@(d2,d3,a1prime,a1,a2,semiz,z,e,w,kappa_j) w*kappa_j*d2*a2*semiz*z*e;
jequaloneDist=zeros([n_a,vfoptions.n_semiz,n_z,vfoptions.n_e],'gpuArray');
jequaloneDist(1,1,2,ceil(n_z/2),ceil(vfoptions.n_e/2))=1;

[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

cases={'experienceassetsemiz runs (positive control)','experienceassetsemiz and experienceasset both set','experienceassetsemiz and experienceassetu both set','experienceassetsemiz without the shock it depends on'};
expected={'','only one of','only one of','needs both z and e'};
cmds={'AutoCorrTransProbs_FHorz','CrossSectionCovarCorr_FHorz (M=2)','AgeConditionalStats_CrossSectionCovarCorr_FHorz (M=2)','LifeCycleProfiles_FHorz_Case1 (M=2)'};
nchecks=0;
for kk=1:length(cases)
    for cmd=1:length(cmds)
        if kk==4 && cmd>1 % the shock check of the other three sits after their single-period statistics (see above)
            continue
        end
        so=simoptions;
        so.ntimeperiods=[1,2];
        n_zc=n_z; z_gridc=z_grid; pi_zc=pi_z;
        if kk==2
            so.experienceasset=1;
        elseif kk==3
            so.experienceassetu=1;
            so.n_u=1; so.u_grid=1; so.pi_u=1; % so that it is the two flags that are refused, not missing u options
        elseif kk==4
            so.n_semiz=0; % no semiz
        end
        errmsg='';
        try
            if cmd==1
                EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,F,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,pi_zc,so);
            elseif cmd==2
                EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,F,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,so);
            elseif cmd==3
                EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,F,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,so);
            else
                LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,F,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,so);
            end
        catch ME
            errmsg=ME.message;
        end
        nchecks=nchecks+1;
        if kk==1
            bad=~isempty(errmsg);
            fprintf('MultiPeriodStats semiz controls: %s, %s, errored, this should be zero: %i \n',cmds{cmd},cases{kk},bad)
        else
            bad=~contains(errmsg,expected{kk});
            fprintf('MultiPeriodStats semiz controls: %s, %s, not refused with the expected error, this should be zero: %i \n',cmds{cmd},cases{kk},bad)
        end
        if bad && ~isempty(errmsg)
            fprintf('    the error was: %s \n',errmsg)
        end
    end
end
fprintf('MultiPeriodStats semiz controls: %i (command, case) pairs checked, of %i \n',nchecks,3*length(cmds)+1)

%% PType: two identical types of the experienceassetsemiz model, per type equal to the solo command
Names_i={'typeA','typeB'};
PTypeDistParamNames={'ptypemass'};
ParamsPT=Params;
ParamsPT.ptypemass=[0.5;0.5];
[~,PolicyPT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsPT,DiscountFactorParamNames,vfoptions);
StationaryDistPT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,PolicyPT,n_d,n_a,n_z,N_j,Names_i,pi_z,ParamsPT,simoptions);
soPT=simoptions; soPT.timehorizons=[1,2]; soPT.ntimeperiods=[1,2];
ACPT=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDistPT,PolicyPT,F,ParamsPT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,soPT);
ACsolo=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,F,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,soPT);
fprintf('MultiPeriodStats semiz PType: human has zero StdDeviation at every age (the agents never move, the comparison would be vacuous), this should be zero: %i \n',all(gather(ACsolo.human.StdDeviation)<1e-12 | isnan(gather(ACsolo.human.StdDeviation))))
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
    fprintf('MultiPeriodStats semiz PType (type %s): %i numbers compared, largest difference at %s \n',Names_i{ii},nnumeric,worst)
    fprintf('MultiPeriodStats semiz PType (type %s): AutoCorrTransProbs_FHorz_PType per type vs solo, this should be zero: %.3e \n',Names_i{ii},maxdiff)
    fprintf('MultiPeriodStats semiz PType (type %s): entries NaN in only one, this should be zero: %i \n',Names_i{ii},nnanbad)
    fprintf('MultiPeriodStats semiz PType (type %s): fields missing or of a different size, this should be zero: %i \n',Names_i{ii},nmismatch)
    fprintf('MultiPeriodStats semiz PType (type %s): nothing compared (a vacuous pass), this should be zero: %i \n',Names_i{ii},nnumeric==0)
end

output=struct();

end
