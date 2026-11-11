function output=CoreFHorzExpAsset_MultiPeriodStats_FakeExpAsset(n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline)
% The multi-period statistics with an experience asset that is really a standard endogenous state (the 'fake experience asset' of
% CrossTests3): aprimeFn=@(d2,a2) d2 with d2_grid=a2_grid, so a2'=d2 is always on the grid. The model is solved twice, as a
% standard two-endogenous-state model (a1',a2' chosen) and as an experience-asset model (d2,a1' chosen), and every output of
%   EvalFnOnAgentDist_AutoCorrTransProbs_FHorz (timehorizons [1,2], ntimeperiods [1,2], a conditional restriction, transprobs;
%       and a second call with the log of a 2-period window, timehorizons 2),
%   EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz (ntimeperiods [1,2], a conditional restriction),
%   EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz (the same),
%   LifeCycleProfiles_FHorz_Case1 (the same)
% must be the same. The two models go through different policy steps in these commands: the standard model through the branch
% that reads a' from the Policy, the experience asset through the branch that builds the a2' corners from aprimeFn. So this is an
% independent check of the experience-asset branch (its a' indexing, its corner weights, its z and e offsets).
% Two cases: no z and no e (the N_semizze==0 path), and z and e (the return depends on z and e through the wage, so the policy
% differs across (z,e) and an ordering mistake in the z/e offsets would show).
% Plus a PType cell: two identical types of the experience-asset model (z and e), EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType
% per type must equal the solo command.
% Every comparison is over all numeric fields (recursively, through the conditional-restriction structures and the TransitionProbs
% cells): the largest difference, the number of entries NaN in one but not the other, and the number of fields or cells missing
% or of a different size in one of the two, each printed as 'this should be zero'.

fprintf('\n========== CoreFHorzExpAsset_MultiPeriodStats_FakeExpAsset ==========\n')

aprimeFn=@(d2,a2) d2; % d2 is just a2prime
d_grid=a_grid(n_a(1)+1:end); % d2 grid is the a2 grid
n_d=n_a(2);

for cc=1:2
    if cc==1
        casestr='noz noe';
        n_zc=0;
        z_gridc=[];
        pi_zc=[];
        vfoptions0=struct();
        simoptions0=struct();
        ReturnFn0=@(a1prime,a2prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_TwoEndo_nod1_noz_noe_nosemiz(a1prime,a2prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension);
        ReturnFn1=@(d2,a1prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz(d2,a1prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension);
        F0=struct(); F1=struct(); R0=struct(); R1=struct();
        F0.assets=@(a1prime,a2prime,a1,a2) a1;
        F0.human=@(a1prime,a2prime,a1,a2) a2;
        F0.earnings=@(a1prime,a2prime,a1,a2,w,kappa_j) w*kappa_j*a2prime*a2;
        F0.assetsplus=@(a1prime,a2prime,a1,a2) a1prime+1; % positive (for the log), and a function of the decision
        F1.assets=@(d2,a1prime,a1,a2) a1;
        F1.human=@(d2,a1prime,a1,a2) a2;
        F1.earnings=@(d2,a1prime,a1,a2,w,kappa_j) w*kappa_j*d2*a2;
        F1.assetsplus=@(d2,a1prime,a1,a2) a1prime+1;
        R0.highhuman=@(a1prime,a2prime,a1,a2) (a2>1); % human capital is chosen freely (a2'=d2) and raises earnings, so many agents have a2>1 (a1>0.5 had zero mass at every age)
        R1.highhuman=@(d2,a1prime,a1,a2) (a2>1); % human capital is chosen freely (a2'=d2) and raises earnings, so many agents have a2>1 (a1>0.5 had zero mass at every age)
        jequaloneDist=zeros(n_a,'gpuArray');
        jequaloneDist(1,2:end)=1/(n_a(2)-1); % no a1, spread over a2>0 (at a2=0 earnings are zero and every choice is -Inf; and with no shocks a single starting point would leave no cross-section at all)
    else
        casestr='z e';
        n_zc=n_z;
        z_gridc=z_grid;
        pi_zc=pi_z;
        vfoptions0=struct();
        vfoptions0.n_e=vfoptionsbaseline.n_e;
        vfoptions0.e_grid=vfoptionsbaseline.e_grid;
        vfoptions0.pi_e=vfoptionsbaseline.pi_e;
        simoptions0=vfoptions0;
        simoptions0.pi_z=pi_z; % the M>1 windows of the cross-section commands and LifeCycleProfiles need the z transitions
        % the wage is w*z*e, so the policy depends on (z,e)
        ReturnFn0=@(a1prime,a2prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_TwoEndo_nod1_noz_noe_nosemiz(a1prime,a2prime,a1,a2,r,w*z*e,kappa_j,sigma,agej,Jr,pension);
        ReturnFn1=@(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz(d2,a1prime,a1,a2,r,w*z*e,kappa_j,sigma,agej,Jr,pension);
        F0=struct(); F1=struct(); R0=struct(); R1=struct();
        F0.assets=@(a1prime,a2prime,a1,a2,z,e) a1;
        F0.human=@(a1prime,a2prime,a1,a2,z,e) a2;
        F0.earnings=@(a1prime,a2prime,a1,a2,z,e,w,kappa_j) w*kappa_j*a2prime*a2*z*e;
        F0.assetsplus=@(a1prime,a2prime,a1,a2,z,e) a1prime+1;
        F1.assets=@(d2,a1prime,a1,a2,z,e) a1;
        F1.human=@(d2,a1prime,a1,a2,z,e) a2;
        F1.earnings=@(d2,a1prime,a1,a2,z,e,w,kappa_j) w*kappa_j*d2*a2*z*e;
        F1.assetsplus=@(d2,a1prime,a1,a2,z,e) a1prime+1;
        R0.highhuman=@(a1prime,a2prime,a1,a2,z,e) (a2>1); % human capital is chosen freely (a2'=d2) and raises earnings, so many agents have a2>1 (a1>0.5 had zero mass at every age)
        R1.highhuman=@(d2,a1prime,a1,a2,z,e) (a2>1); % human capital is chosen freely (a2'=d2) and raises earnings, so many agents have a2>1 (a1>0.5 had zero mass at every age)
        jequaloneDist=zeros([n_a,n_z,vfoptions0.n_e],'gpuArray');
        jequaloneDist(1,ceil(n_a(2)/2),ceil(n_z/2),ceil(vfoptions0.n_e/2))=1; % no a1, a2 mid-grid (at a2=0 every choice is -Inf)
    end
    vfoptions1=vfoptions0;
    vfoptions1.experienceasset=1;
    vfoptions1.aprimeFn=aprimeFn;
    simoptions1=simoptions0;
    simoptions1.experienceasset=1;
    simoptions1.aprimeFn=aprimeFn;
    simoptions1.d_grid=d_grid;
    simoptions1.a_grid=a_grid;

    [V0,Policy0]=ValueFnIter_Case1_FHorz(0,n_a,n_zc,N_j,[],a_grid,z_gridc,pi_zc,ReturnFn0,Params,DiscountFactorParamNames,[],vfoptions0);
    StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy0,0,n_a,n_zc,N_j,pi_zc,Params,simoptions0);
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,pi_zc,ReturnFn1,Params,DiscountFactorParamNames,[],vfoptions1);
    StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1,n_d,n_a,n_zc,N_j,pi_zc,Params,simoptions1);
    Policy0alt=[Policy0(2,:); Policy0(1,:)]; % (a2',a1') is (d2,a1')
    fprintf('MultiPeriodStats fake expasset (%s): the two models, V, this should be zero: %.3e \n',casestr,max(abs(V0(:)-V1(:))))
    fprintf('MultiPeriodStats fake expasset (%s): the two models, Policy, this should be zero: %.3e \n',casestr,max(abs(Policy0alt(:)-Policy1(:))))
    fprintf('MultiPeriodStats fake expasset (%s): the two models, StationaryDist, this should be zero: %.3e \n',casestr,max(abs(StationaryDist0(:)-StationaryDist1(:))))

    % AutoCorrTransProbs: horizons 1 and 2, windows of 1 and 2 periods, a restriction, the transition probabilities of human
    soA0=simoptions0; soA0.timehorizons=[1,2]; soA0.ntimeperiods=[1,2]; soA0.conditionalrestrictions=R0; soA0.transprobs={'human'};
    soA1=simoptions1; soA1.timehorizons=[1,2]; soA1.ntimeperiods=[1,2]; soA1.conditionalrestrictions=R1; soA1.transprobs={'human'};
    AC0=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist0,Policy0,F0,Params,[],0,n_a,n_zc,N_j,[],a_grid,z_gridc,pi_zc,soA0);
    AC1=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist1,Policy1,F1,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,pi_zc,soA1);
    fprintf('MultiPeriodStats fake expasset (%s): the restriction highhuman has zero mass at every age (the restricted comparisons would be vacuous), this should be zero: %i \n',casestr,all(gather(AC1.highhuman.RestrictedSampleMass)==0))
    fprintf('MultiPeriodStats fake expasset (%s): assets and human have zero StdDeviation at every age (the agents never move, the comparisons would be vacuous), this should be zero: %i \n',casestr,all(gather(AC1.assets.StdDeviation)<1e-12 | isnan(gather(AC1.assets.StdDeviation))) && all(gather(AC1.human.StdDeviation)<1e-12 | isnan(gather(AC1.human.StdDeviation))))
    % AutoCorrTransProbs: the log of the 2-period window of assetsplus (puts the window on the edges of the transition), horizon 2
    logwc=zeros(4,N_j,2); logwc(4,:,:)=1; % row 4 is assetsplus; unrestricted and restricted
    soL0=simoptions0; soL0.timehorizons=2; soL0.ntimeperiods=2; soL0.conditionalrestrictions=R0; soL0.logwhichcombos=logwc;
    soL1=simoptions1; soL1.timehorizons=2; soL1.ntimeperiods=2; soL1.conditionalrestrictions=R1; soL1.logwhichcombos=logwc;
    ACL0=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist0,Policy0,F0,Params,[],0,n_a,n_zc,N_j,[],a_grid,z_gridc,pi_zc,soL0);
    ACL1=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist1,Policy1,F1,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,pi_zc,soL1);
    % The cross-section commands and LifeCycleProfiles: windows of 1 and 2 periods, a restriction
    soC0=simoptions0; soC0.ntimeperiods=[1,2]; soC0.conditionalrestrictions=R0;
    soC1=simoptions1; soC1.ntimeperiods=[1,2]; soC1.conditionalrestrictions=R1;
    CS0=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist0,Policy0,F0,Params,[],0,n_a,n_zc,N_j,[],a_grid,z_gridc,soC0);
    CS1=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist1,Policy1,F1,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,soC1);
    ACS0=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist0,Policy0,F0,Params,[],0,n_a,n_zc,N_j,[],a_grid,z_gridc,soC0);
    ACS1=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist1,Policy1,F1,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,soC1);
    LCP0=LifeCycleProfiles_FHorz_Case1(StationaryDist0,Policy0,F0,Params,[],0,n_a,n_zc,N_j,[],a_grid,z_gridc,soC0);
    LCP1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1,F1,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,soC1);

    pairs={AC0,AC1,'AutoCorrTransProbs'; ACL0,ACL1,'AutoCorrTransProbs (log of 2-period window)'; CS0,CS1,'CrossSectionCovarCorr'; ACS0,ACS1,'AgeConditionalStats_CrossSectionCovarCorr'; LCP0,LCP1,'LifeCycleProfiles'};
    for pp=1:size(pairs,1)
        % Walk the two outputs together: every numeric field, through the sub-structures and cells
        stack={{pairs{pp,1},pairs{pp,2},''}};
        maxdiff=0; worst='(none)'; nnanbad=0; nmismatch=0; nnumeric=0;
        while ~isempty(stack)
            item=stack{end}; stack(end)=[];
            X=item{1}; Y=item{2};
            if isstruct(X)
                if ~isstruct(Y)
                    nmismatch=nmismatch+1;
                    continue
                end
                fnX=fieldnames(X); fnY=fieldnames(Y);
                nmismatch=nmismatch+length(setxor(fnX,fnY));
                fnXY=intersect(fnX,fnY);
                for kk=1:length(fnXY)
                    stack{end+1}={X.(fnXY{kk}),Y.(fnXY{kk}),[item{3},'.',fnXY{kk}]}; %#ok<AGROW>
                end
            elseif iscell(X)
                if ~iscell(Y) || ~isequal(size(X),size(Y))
                    nmismatch=nmismatch+1;
                    continue
                end
                for kk=1:numel(X)
                    stack{end+1}={X{kk},Y{kk},[item{3},'{',num2str(kk),'}']}; %#ok<AGROW>
                end
            elseif isnumeric(X) || islogical(X)
                Xd=gather(double(X)); Yd=gather(double(Y));
                if ~isequal(size(Xd),size(Yd))
                    nmismatch=nmismatch+1;
                    continue
                end
                nnumeric=nnumeric+numel(Xd);
                % NaN in one and a number in the other. The commands report a correlation only where the product of the standard
                % deviations exceeds 1e-15; where the (restricted) population has no variance that product is double noise, which
                % falls on either side of the cutoff depending on the order of the sums, and the 'correlation' reported is itself
                % noise (~1e-16). So a NaN against a number below 1e-12 is the same answer (zero/undefined), printed but not
                % counted; a NaN against anything larger is a real difference.
                nanbad=find(isnan(Xd(:))~=isnan(Yd(:)));
                for kk=1:length(nanbad) % say where, so a NaN-pattern difference can be read off the diary
                    defined=[Xd(nanbad(kk)),Yd(nanbad(kk))];
                    defined=defined(~isnan(defined));
                    if abs(defined)>1e-12
                        nnanbad=nnanbad+1;
                        fprintf('    NaN in only one (counted): %s, entry %i: standard asset %.6e, experience asset %.6e \n',item{3},nanbad(kk),Xd(nanbad(kk)),Yd(nanbad(kk)))
                    else
                        fprintf('    NaN in only one, against double noise (not counted): %s, entry %i: standard asset %.6e, experience asset %.6e \n',item{3},nanbad(kk),Xd(nanbad(kk)),Yd(nanbad(kk)))
                    end
                end
                both=~isnan(Xd(:)) & ~isnan(Yd(:));
                if any(both)
                    dd=max(abs(Xd(both)-Yd(both)));
                    if dd>maxdiff
                        maxdiff=dd;
                        worst=item{3};
                    end
                end
            end % char (the Notes) is not compared
        end
        fprintf('MultiPeriodStats fake expasset (%s): %s, %i numbers compared, largest difference at %s \n',casestr,pairs{pp,3},nnumeric,worst)
        fprintf('MultiPeriodStats fake expasset (%s): %s vs standard asset, this should be zero: %.3e \n',casestr,pairs{pp,3},maxdiff)
        fprintf('MultiPeriodStats fake expasset (%s): %s, entries NaN in only one, this should be zero: %i \n',casestr,pairs{pp,3},nnanbad)
        fprintf('MultiPeriodStats fake expasset (%s): %s, fields missing or of a different size, this should be zero: %i \n',casestr,pairs{pp,3},nmismatch)
        fprintf('MultiPeriodStats fake expasset (%s): %s, nothing compared (a vacuous pass), this should be zero: %i \n',casestr,pairs{pp,3},nnumeric==0)
    end

    if cc==2
        %% PType: two identical types of the experience-asset model, per type equal to the solo command
        Names_i={'typeA','typeB'};
        PTypeDistParamNames={'ptypemass'};
        ParamsPT=Params;
        ParamsPT.ptypemass=[0.5;0.5];
        [~,PolicyPT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_zc,N_j,Names_i,d_grid,a_grid,z_gridc,pi_zc,ReturnFn1,ParamsPT,DiscountFactorParamNames,vfoptions1);
        StationaryDistPT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,PolicyPT,n_d,n_a,n_zc,N_j,Names_i,pi_zc,ParamsPT,simoptions1);
        soPT=simoptions1; soPT.timehorizons=[1,2]; soPT.ntimeperiods=[1,2];
        ACPT=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDistPT,PolicyPT,F1,ParamsPT,n_d,n_a,n_zc,N_j,Names_i,d_grid,a_grid,z_gridc,pi_zc,soPT);
        ACsolo=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist1,Policy1,F1,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,pi_zc,soPT);
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
                    nanbad=find(isnan(Xd(:))~=isnan(Yd(:))); % a NaN against double noise (below 1e-12) is the same answer, as above
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
            fprintf('MultiPeriodStats fake expasset PType (type %s): %i numbers compared, largest difference at %s \n',Names_i{ii},nnumeric,worst)
            fprintf('MultiPeriodStats fake expasset PType (type %s): AutoCorrTransProbs_FHorz_PType per type vs solo, this should be zero: %.3e \n',Names_i{ii},maxdiff)
            fprintf('MultiPeriodStats fake expasset PType (type %s): entries NaN in only one, this should be zero: %i \n',Names_i{ii},nnanbad)
            fprintf('MultiPeriodStats fake expasset PType (type %s): fields missing or of a different size, this should be zero: %i \n',Names_i{ii},nmismatch)
            fprintf('MultiPeriodStats fake expasset PType (type %s): nothing compared (a vacuous pass), this should be zero: %i \n',Names_i{ii},nnumeric==0)
        end
    end
end

output=struct();

end
