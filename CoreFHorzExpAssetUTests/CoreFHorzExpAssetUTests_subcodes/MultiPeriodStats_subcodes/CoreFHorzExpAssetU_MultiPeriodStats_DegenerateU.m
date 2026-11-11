function output=CoreFHorzExpAssetU_MultiPeriodStats_DegenerateU(n_d2,d2_grid,n_a,a_grid,n_z,N_j,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline)
% The multi-period statistics with experienceassetu and a u that is really no shock: n_u=2, u_grid=[1;1], pi_u=[0.3;0.7]. With the
% baseline aprimeFn u*(phi1*(1-d2)+(1-phi2)*a2) that is the experienceasset with aprimeFn phi1*(1-d2)+(1-phi2)*a2, so the model is
% solved both ways (experienceassetu, experienceasset) and every output of
%   EvalFnOnAgentDist_AutoCorrTransProbs_FHorz (timehorizons [1,2], ntimeperiods [1,2], a conditional restriction, transprobs;
%       and a second call with the log of a 2-period window, timehorizons 2),
%   EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz (ntimeperiods [1,2], a conditional restriction),
%   EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz (the same),
%   LifeCycleProfiles_FHorz_Case1 (the same)
% must be the same. The two u values are equal but their probabilities are not, so the experienceassetu side carries every
% a2' corner twice, with weights 0.3 and 0.7: a mistake in the u columns of the policy step (their order, their pi_u weights, the
% offsets of the corners within them) would show. A genuine u is checked against a brute force in
% CoreFHorzExpAssetU_MultiPeriodStats_BruteForce.
% Two cases: no z and no e (the N_semizze==0 path) and z and e (the return depends on z and e through the wage).
% Plus a PType cell: two identical types of the experienceassetu model (z and e), EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType
% per type must equal the solo command.
% Every comparison is over all numeric fields (recursively): the largest difference, the number of entries NaN in one but not the
% other (against a number above 1e-12: the commands report a correlation only where the product of the standard deviations
% exceeds 1e-15, and where a population has no variance that product is double noise on either side of the cutoff), and the
% number of fields or cells missing or of a different size, each printed as 'this should be zero'.

fprintf('\n========== CoreFHorzExpAssetU_MultiPeriodStats_DegenerateU ==========\n')

aprimeFn_expasset=@(d2,a2,phi1,phi2) phi1*(1-d2)+(1-phi2)*a2; % = the baseline experienceassetu aprimeFn at u=1
n_d=n_d2;
d_grid=d2_grid;

for cc=1:2
    if cc==1
        casestr='noz noe';
        n_zc=0;
        z_gridc=[];
        pi_zc=[];
        vfoptions0=struct();
        simoptions0=struct();
        ReturnFn=@(d2,a1prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz(d2,a1prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension);
        F=struct(); R=struct();
        F.assets=@(d2,a1prime,a1,a2) a1;
        F.human=@(d2,a1prime,a1,a2) a2;
        F.earnings=@(d2,a1prime,a1,a2,w,kappa_j) w*kappa_j*d2*a2;
        F.assetsplus=@(d2,a1prime,a1,a2) a1prime+1; % positive (for the log), and a function of the decision
        R.highhuman=@(d2,a1prime,a1,a2) (a2>1);
        jequaloneDist=zeros(n_a,'gpuArray');
        jequaloneDist(1,2:end)=1/(n_a(2)-1); % no a1, spread over a2>0 (with no shocks a single starting point would leave no cross-section)
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
        ReturnFn=@(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_e_nosemiz(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
        F=struct(); R=struct();
        F.assets=@(d2,a1prime,a1,a2,z,e) a1;
        F.human=@(d2,a1prime,a1,a2,z,e) a2;
        F.earnings=@(d2,a1prime,a1,a2,z,e,w,kappa_j) w*kappa_j*d2*a2*z*e;
        F.assetsplus=@(d2,a1prime,a1,a2,z,e) a1prime+1;
        R.highhuman=@(d2,a1prime,a1,a2,z,e) (a2>1);
        jequaloneDist=zeros([n_a,n_z,vfoptions0.n_e],'gpuArray');
        jequaloneDist(1,ceil(n_a(2)/2),ceil(n_z/2),ceil(vfoptions0.n_e/2))=1; % no a1, a2 mid-grid
    end
    % Model 0: experienceasset
    vfoptionsA=vfoptions0;
    vfoptionsA.experienceasset=1;
    vfoptionsA.aprimeFn=aprimeFn_expasset;
    simoptionsA=simoptions0;
    simoptionsA.experienceasset=1;
    simoptionsA.aprimeFn=aprimeFn_expasset;
    simoptionsA.d_grid=d_grid;
    simoptionsA.a_grid=a_grid;
    % Model 1: experienceassetu, two equal u values with unequal probabilities
    vfoptionsB=vfoptions0;
    vfoptionsB.experienceassetu=1;
    vfoptionsB.aprimeFn=vfoptionsbaseline.aprimeFn; % u*(phi1*(1-d2)+(1-phi2)*a2)
    vfoptionsB.n_u=2;
    vfoptionsB.u_grid=[1;1];
    vfoptionsB.pi_u=[0.3;0.7];
    simoptionsB=simoptions0;
    simoptionsB.experienceassetu=1;
    simoptionsB.aprimeFn=vfoptionsB.aprimeFn;
    simoptionsB.n_u=vfoptionsB.n_u;
    simoptionsB.u_grid=vfoptionsB.u_grid;
    simoptionsB.pi_u=vfoptionsB.pi_u;
    simoptionsB.d_grid=d_grid;
    simoptionsB.a_grid=a_grid;

    [V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,pi_zc,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsA);
    StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy0,n_d,n_a,n_zc,N_j,pi_zc,Params,simoptionsA);
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,pi_zc,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsB);
    StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1,n_d,n_a,n_zc,N_j,pi_zc,Params,simoptionsB);
    fprintf('MultiPeriodStats degenerate u (%s): the two models, V, this should be zero: %.3e \n',casestr,max(abs(V0(:)-V1(:))))
    fprintf('MultiPeriodStats degenerate u (%s): the two models, Policy, this should be zero: %.3e \n',casestr,max(abs(Policy0(:)-Policy1(:))))
    fprintf('MultiPeriodStats degenerate u (%s): the two models, StationaryDist, this should be zero: %.3e \n',casestr,max(abs(StationaryDist0(:)-StationaryDist1(:))))

    % AutoCorrTransProbs: horizons 1 and 2, windows of 1 and 2 periods, a restriction, the transition probabilities of human
    soA0=simoptionsA; soA0.timehorizons=[1,2]; soA0.ntimeperiods=[1,2]; soA0.conditionalrestrictions=R; soA0.transprobs={'human'};
    soA1=simoptionsB; soA1.timehorizons=[1,2]; soA1.ntimeperiods=[1,2]; soA1.conditionalrestrictions=R; soA1.transprobs={'human'};
    AC0=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist0,Policy0,F,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,pi_zc,soA0);
    AC1=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist1,Policy1,F,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,pi_zc,soA1);
    fprintf('MultiPeriodStats degenerate u (%s): the restriction highhuman has zero mass at every age (the restricted comparisons would be vacuous), this should be zero: %i \n',casestr,all(gather(AC1.highhuman.RestrictedSampleMass)==0))
    fprintf('MultiPeriodStats degenerate u (%s): assets and human have zero StdDeviation at every age (the agents never move, the comparisons would be vacuous), this should be zero: %i \n',casestr,all(gather(AC1.assets.StdDeviation)<1e-12 | isnan(gather(AC1.assets.StdDeviation))) && all(gather(AC1.human.StdDeviation)<1e-12 | isnan(gather(AC1.human.StdDeviation))))
    % AutoCorrTransProbs: the log of the 2-period window of assetsplus (puts the window on the edges of the transition), horizon 2
    logwc=zeros(4,N_j,2); logwc(4,:,:)=1; % row 4 is assetsplus; unrestricted and restricted
    soL0=simoptionsA; soL0.timehorizons=2; soL0.ntimeperiods=2; soL0.conditionalrestrictions=R; soL0.logwhichcombos=logwc;
    soL1=simoptionsB; soL1.timehorizons=2; soL1.ntimeperiods=2; soL1.conditionalrestrictions=R; soL1.logwhichcombos=logwc;
    ACL0=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist0,Policy0,F,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,pi_zc,soL0);
    ACL1=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist1,Policy1,F,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,pi_zc,soL1);
    % The cross-section commands and LifeCycleProfiles: windows of 1 and 2 periods, a restriction
    soC0=simoptionsA; soC0.ntimeperiods=[1,2]; soC0.conditionalrestrictions=R;
    soC1=simoptionsB; soC1.ntimeperiods=[1,2]; soC1.conditionalrestrictions=R;
    CS0=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist0,Policy0,F,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,soC0);
    CS1=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist1,Policy1,F,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,soC1);
    ACS0=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist0,Policy0,F,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,soC0);
    ACS1=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist1,Policy1,F,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,soC1);
    LCP0=LifeCycleProfiles_FHorz_Case1(StationaryDist0,Policy0,F,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,soC0);
    LCP1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1,F,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,soC1);

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
                nanbad=find(isnan(Xd(:))~=isnan(Yd(:)));
                for kk=1:length(nanbad) % say where, so a NaN-pattern difference can be read off the diary
                    defined=[Xd(nanbad(kk)),Yd(nanbad(kk))];
                    defined=defined(~isnan(defined));
                    if abs(defined)>1e-12
                        nnanbad=nnanbad+1;
                        fprintf('    NaN in only one (counted): %s, entry %i: experienceasset %.6e, experienceassetu %.6e \n',item{3},nanbad(kk),Xd(nanbad(kk)),Yd(nanbad(kk)))
                    else
                        fprintf('    NaN in only one, against double noise (not counted): %s, entry %i: experienceasset %.6e, experienceassetu %.6e \n',item{3},nanbad(kk),Xd(nanbad(kk)),Yd(nanbad(kk)))
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
        fprintf('MultiPeriodStats degenerate u (%s): %s, %i numbers compared, largest difference at %s \n',casestr,pairs{pp,3},nnumeric,worst)
        fprintf('MultiPeriodStats degenerate u (%s): %s, experienceassetu vs experienceasset, this should be zero: %.3e \n',casestr,pairs{pp,3},maxdiff)
        fprintf('MultiPeriodStats degenerate u (%s): %s, entries NaN in only one, this should be zero: %i \n',casestr,pairs{pp,3},nnanbad)
        fprintf('MultiPeriodStats degenerate u (%s): %s, fields missing or of a different size, this should be zero: %i \n',casestr,pairs{pp,3},nmismatch)
        fprintf('MultiPeriodStats degenerate u (%s): %s, nothing compared (a vacuous pass), this should be zero: %i \n',casestr,pairs{pp,3},nnumeric==0)
    end

    if cc==2
        %% PType: two identical types of the experienceassetu model, per type equal to the solo command
        Names_i={'typeA','typeB'};
        PTypeDistParamNames={'ptypemass'};
        ParamsPT=Params;
        ParamsPT.ptypemass=[0.5;0.5];
        [~,PolicyPT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_zc,N_j,Names_i,d_grid,a_grid,z_gridc,pi_zc,ReturnFn,ParamsPT,DiscountFactorParamNames,vfoptionsB);
        StationaryDistPT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,PolicyPT,n_d,n_a,n_zc,N_j,Names_i,pi_zc,ParamsPT,simoptionsB);
        soPT=simoptionsB; soPT.timehorizons=[1,2]; soPT.ntimeperiods=[1,2];
        ACPT=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDistPT,PolicyPT,F,ParamsPT,n_d,n_a,n_zc,N_j,Names_i,d_grid,a_grid,z_gridc,pi_zc,soPT);
        ACsolo=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist1,Policy1,F,Params,[],n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,pi_zc,soPT);
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
            fprintf('MultiPeriodStats degenerate u PType (type %s): %i numbers compared, largest difference at %s \n',Names_i{ii},nnumeric,worst)
            fprintf('MultiPeriodStats degenerate u PType (type %s): AutoCorrTransProbs_FHorz_PType per type vs solo, this should be zero: %.3e \n',Names_i{ii},maxdiff)
            fprintf('MultiPeriodStats degenerate u PType (type %s): entries NaN in only one, this should be zero: %i \n',Names_i{ii},nnanbad)
            fprintf('MultiPeriodStats degenerate u PType (type %s): fields missing or of a different size, this should be zero: %i \n',Names_i{ii},nmismatch)
            fprintf('MultiPeriodStats degenerate u PType (type %s): nothing compared (a vacuous pass), this should be zero: %i \n',Names_i{ii},nnumeric==0)
        end
    end
end

output=struct();

end
