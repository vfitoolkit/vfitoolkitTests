function output=CoreFHorzTPathPType_PerTypeFnsToEvaluate(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,Names_i,transpathoptionsbaseline)
% Test 5: two PTypes that are actually identical. Solve once, then evaluate
% FnsToEvaluate two ways: (A) a single function returning the endogenous
% state a, broadcast to both types; (B) per-type functions (struct-of-structs,
% one sub-entry per type under the same field name, together covering all
% agents) each returning a. Both must give the same answer.
% Checked for EvalFnOnTransPath_AggVars_Case1_FHorz_PType (aggregate and
% per-type breakdowns), for one iteration of TransitionPath_Case1_FHorz_PType
% (shooting, and Anderson), and for EvalFnOnTransPath_AllStats_Case1_FHorz_PType
% (also with conditionalrestrictions, which is checked against the stationary
% PType AllStats period by period) and LifeCycleProfiles_TransPath_FHorz_Case1_PType.
% Only z is used (no d, no e, no semiz).

n_d=0;
d_grid=[];

% ReturnFn declares kappa_j (not kappa_j_pt), so both types are identical.
ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

% (A) single FnsToEvaluate, broadcast to both types
FnsToEvaluate_single.assets=@(aprime,a,z) a;
% (B) per-type FnsToEvaluate: same field name, one sub-entry per type
FnsToEvaluate_pertype.assets.(Names_i{1})=@(aprime,a,z) a;
FnsToEvaluate_pertype.assets.(Names_i{2})=@(aprime,a,z) a;

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

vfoptions=struct();
simoptions=struct();

% Trivial period-T V and Policy (as in CoreFHorzTPathTests)
V_final_PT=struct(); Policy_final_PT=struct();
for ii=1:N_i
    V_final_PT.(Names_i{ii})=zeros([n_a,n_z,N_j],'gpuArray');
    Policy_final_PT.(Names_i{ii})=ones([1,n_a,n_z,N_j],'gpuArray');
end

%% Solve the path machinery once (FnsToEvaluate does not enter these, so both cases share it)
[~,PolicyPath_PT]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, Policy_final_PT, Params, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
Policy_init_PT=struct();
for ii=1:N_i
    Policy_init_PT.(Names_i{ii})=PolicyPath_PT.(Names_i{ii})(:,:,:,:,1);
end
AgentDist_initial_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_init_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AgentDistPath_PT=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist, PricePath, ParamPath, PolicyPath_PT, AgeWeightParamNames, PTypeDistParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z, T,Params, transpathoptionsbaseline, simoptions);

%% AggVars on the transition path: single vs per-type FnsToEvaluate
AggVarsPath_A=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate_single, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
AggVarsPath_B=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate_pertype, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);

fprintf('per-type FnsToEvaluate, AggVarsPath assets.Mean, this should be zero: %.3e \n',max(abs(AggVarsPath_A.assets.Mean(:)-AggVarsPath_B.assets.Mean(:))))
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('per-type FnsToEvaluate, AggVarsPath assets breakdown (type %s), this should be zero: %.3e \n',nm,max(abs(AggVarsPath_A.assets.(nm).Mean(:)-AggVarsPath_B.assets.(nm).Mean(:))))
end

%% One iteration of the GE transition path: single vs per-type FnsToEvaluate
% maxiter=2: since the shooting algorithm returns the path whose general eqm conditions it evaluated,
% maxiter=1 would return the starting path unchanged and compare it with itself. The general eqm
% condition depends on assets, so the update does too.
transpathoptions=transpathoptionsbaseline;
transpathoptions.maxiter=2;
transpathoptions.GEnewprice=3; % the shooting algorithm
transpathoptions.GEnewprice3.howtoupdate={'dummy','r',0,0.1};

GeneralEqmEqns.dummy=@(r,assets) r-0.02-0.01*assets;

[PricePath_A,GEcondnPath_A]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_single, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
[PricePath_B,GEcondnPath_B]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_pertype, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);

fprintf('per-type FnsToEvaluate, one iter of GE TPath, this should be zero: %.3e \n',max(abs(PricePath_A.r(:)-PricePath_B.r(:))))
fprintf('per-type FnsToEvaluate, one iter of GE TPath, GEcondnPath, this should be zero: %.3e \n',max(abs(GEcondnPath_A.dummy(:)-GEcondnPath_B.dummy(:))))
fprintf('per-type FnsToEvaluate, one iter of GE TPath, for scale the update moved r by %.3e (must not be zero, or the two checks above are empty) \n',max(abs(PricePath_A.r(:)-PricePath.r(:))))

%% One iteration of the GE transition path again, with Anderson acceleration (GEnewprice=2)
% warmup=0 and maxiter=3, so iterations 2 and 3 are Anderson steps (at maxiter=2 Anderson would still be
% the shooting algorithm here, as iteration 1 has no history to mix)
transpathoptionsA=rmfield(transpathoptions,'GEnewprice3');
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=transpathoptions.GEnewprice3.howtoupdate;
transpathoptionsA.maxiter=3;
transpathoptionsA.anderson.warmup=0;
transpathoptionsA.anderson.safeguard=0;
[PricePathA_A,GEcondnPathA_A]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_single, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsA, simoptions, vfoptions);
[PricePathA_B,GEcondnPathA_B]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_pertype, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsA, simoptions, vfoptions);
transpathoptionsS3=transpathoptions;
transpathoptionsS3.maxiter=3;
PricePathS3_A=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_single, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsS3, simoptions, vfoptions);

fprintf('per-type FnsToEvaluate, Anderson at maxiter=3, GE TPath, this should be zero: %.3e \n',max(abs(PricePathA_A.r(:)-PricePathA_B.r(:))))
fprintf('per-type FnsToEvaluate, Anderson at maxiter=3, GE TPath, GEcondnPath, this should be zero: %.3e \n',max(abs(GEcondnPathA_A.dummy(:)-GEcondnPathA_B.dummy(:))))
fprintf('per-type FnsToEvaluate, Anderson vs shooting, both at maxiter=3, this should NOT be zero: %.3e \n',max(abs(PricePathA_A.r(:)-PricePathS3_A.r(:))))

%% AllStats and life-cycle profiles on the transition path: single vs per-type FnsToEvaluate
% Every field of the single-FnsToEvaluate output (pooled, by ptype, to any depth) must be in the
% per-type output with the same value. With simoptions.conditionalrestrictions too, and the restricted
% AllStats is also checked against the stationary PType AllStats period by period.
OutputTypes={'Mean','Median','Variance','LorenzCurve','Gini','QuantileCutoffs','QuantileMeans'};
simoptionsR=simoptions;
simoptionsR.conditionalrestrictions.wealthy=@(aprime,a,z) (a>1);
AllStatsPath_A=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsToEvaluate_single, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
AllStatsPath_B=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsToEvaluate_pertype, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
AllStatsPathR_A=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsToEvaluate_single, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptionsR);
AllStatsPathR_B=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsToEvaluate_pertype, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptionsR);
cmpcases={'AllStatsPath', AllStatsPath_A, AllStatsPath_B;
          'AllStatsPath with conditionalrestrictions', AllStatsPathR_A, AllStatsPathR_B};
for cc=1:size(cmpcases,1)
    REF=cmpcases{cc,2};
    OUT=cmpcases{cc,3};
    maxdiff=0; nanmismatch=0;
    queue={{}};
    while ~isempty(queue)
        fieldpath=queue{1};
        queue(1)=[];
        if isempty(fieldpath)
            O=REF;
        else
            O=getfield(REF,fieldpath{:});
        end
        if isstruct(O)
            Onames=fieldnames(O);
            for kk=1:length(Onames)
                queue{end+1}=[fieldpath,Onames(kk)];
            end
        else
            P=OUT;
            found=1;
            for kk=1:length(fieldpath)
                if isstruct(P) && isfield(P,fieldpath{kk})
                    P=P.(fieldpath{kk});
                else
                    found=0;
                    break
                end
            end
            if found==0
                nanmismatch=nanmismatch+1;
            elseif iscell(O) % text (the comments): must be the same text
                nanmismatch=nanmismatch+(~isequal(P,O));
            else
                maxdiff=max([maxdiff; abs(P(:)-O(:))],[],'omitnan');
                nanmismatch=nanmismatch+sum(isnan(P(:))~=isnan(O(:)));
            end
        end
    end
    fprintf('per-type FnsToEvaluate, %s, every stat (pooled and by ptype), this should be zero: %.3e \n',cmpcases{cc,1},maxdiff)
    fprintf('per-type FnsToEvaluate, %s, every stat (pooled and by ptype), positions NaN (or text) in one but not the other, or missing, this should be zero: %i \n',cmpcases{cc,1},nanmismatch)
end

% The restricted AllStats against the stationary PType AllStats, period by period
PricePathNamesS=fieldnames(PricePath);
ParamPathNamesS=fieldnames(ParamPath);
maxdiff=0; nanmismatch=0; nfields=0;
for tt=1:T
    Params_tt=Params;
    for kk=1:length(PricePathNamesS)
        Params_tt.(PricePathNamesS{kk})=PricePath.(PricePathNamesS{kk})(tt);
    end
    for kk=1:length(ParamPathNamesS)
        Params_tt.(ParamPathNamesS{kk})=ParamPath.(ParamPathNamesS{kk})(tt);
    end
    StationaryDist_tt=struct(); Policy_tt=struct();
    for ii=1:N_i
        StationaryDist_tt.(Names_i{ii})=AgentDistPath_PT.(Names_i{ii})(:,:,:,tt);
        Policy_tt.(Names_i{ii})=PolicyPath_PT.(Names_i{ii})(:,:,:,:,tt);
    end
    StationaryDist_tt.ptweights=AgentDistPath_PT.ptweights;
    AS=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_tt,Policy_tt,FnsToEvaluate_single,Params_tt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptionsR);
    queue={{}};
    while ~isempty(queue)
        fieldpath=queue{1};
        queue(1)=[];
        if isempty(fieldpath)
            O=AS;
        else
            O=getfield(AS,fieldpath{:});
        end
        if isstruct(O)
            Onames=fieldnames(O);
            for kk=1:length(Onames)
                queue{end+1}=[fieldpath,Onames(kk)];
            end
        elseif ischar(O) || isstring(O) || iscell(O) % text (the comments): must be the same text
            P=getfield(AllStatsPathR_A,fieldpath{:});
            nanmismatch=nanmismatch+(~isequal(P{tt},O)); nfields=nfields+1;
        else
            O=gather(O);
            P=getfield(AllStatsPathR_A,fieldpath{:});
            P=P(:,tt);
            maxdiff=max([maxdiff; abs(P(:)-O(:))],[],'omitnan');
            nanmismatch=nanmismatch+sum(isnan(P(:))~=isnan(O(:))); nfields=nfields+1;
        end
    end
end
fprintf('per-type FnsToEvaluate, AllStatsPath with conditionalrestrictions vs the stationary PType AllStats period by period, %i stat-periods compared, this should be zero: %.3e \n',nfields,maxdiff)
fprintf('per-type FnsToEvaluate, AllStatsPath with conditionalrestrictions vs the stationary PType AllStats period by period, positions NaN (or text) in one but not the other, this should be zero: %i \n',nanmismatch)
fprintf('per-type FnsToEvaluate, for scale, the restricted sample mass (pooled) moves between %.3f and %.3f over the path \n',min(AllStatsPathR_A.wealthy.RestrictedSampleMass.TotalAllPTypes),max(AllStatsPathR_A.wealthy.RestrictedSampleMass.TotalAllPTypes))

% Life-cycle profiles: single vs per-type FnsToEvaluate, pooled and by ptype, all seven outputs
AgeCondStatsPath_A=LifeCycleProfiles_TransPath_FHorz_Case1_PType(FnsToEvaluate_single, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
AgeCondStatsPath_B=LifeCycleProfiles_TransPath_FHorz_Case1_PType(FnsToEvaluate_pertype, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
groupnames=[{''},Names_i];
maxdiff=0; nanmismatch=0;
for gg=1:length(groupnames)
    if isempty(groupnames{gg})
        SA=AgeCondStatsPath_A.assets;
        SB=AgeCondStatsPath_B.assets;
    else
        SA=AgeCondStatsPath_A.assets.(groupnames{gg});
        SB=AgeCondStatsPath_B.assets.(groupnames{gg});
    end
    for oo=1:length(OutputTypes)
        PA=[SA.initialages.(OutputTypes{oo})(:); SA.bornduringtranstion.(OutputTypes{oo})(:)];
        PB=[SB.initialages.(OutputTypes{oo})(:); SB.bornduringtranstion.(OutputTypes{oo})(:)];
        maxdiff=max([maxdiff; abs(PA-PB)],[],'omitnan');
        nanmismatch=nanmismatch+sum(isnan(PA)~=isnan(PB));
    end
end
fprintf('per-type FnsToEvaluate, AgeCondStatsPath, all seven outputs (pooled and by ptype), this should be zero: %.3e \n',maxdiff)
fprintf('per-type FnsToEvaluate, AgeCondStatsPath, entries NaN in one but not the other, this should be zero: %i \n',nanmismatch)

output=struct();

end
