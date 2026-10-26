function output=CoreFHorzTPathPType_NiTwoTypes_vsTwoSolves(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,transpathoptionsbaseline,figure_c)
% Test 3: PType TPath with N_i=2 genuinely different types should equal
% solving the same transition twice without PType, once per type, and
% stacking. Per-type VPath/PolicyPath/AgentDistPath must match the
% corresponding single-type solves; the aggregated AggVarsPath.Mean must
% equal the ptypeweight-weighted combination of the two per-type means.
% Only z is used (no d, no e, no semiz).

N_i=2;

n_d=0;
d_grid=[];

% Per-type kappa_j (already set up as Params.kappa_j_pt = [2 x N_j])
% PType solve uses kappa_j_pt; single-type solves overwrite Params.kappa_j per type.
ReturnFn_PT=@(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension);
ReturnFn_NoPT=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate_PT.assets=@(aprime,a,z) a;
FnsToEvaluate_PT.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate_NoPT.assets=@(aprime,a,z) a;
FnsToEvaluate_NoPT.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

vfoptions=struct();
simoptions=struct();

% Trivial period-T V and Policy (as in CoreFHorzTPathTests)
V_final=zeros([n_a,n_z,N_j],'gpuArray');
Policy_final=ones([1,n_a,n_z,N_j],'gpuArray');
names_PT={'ptype001','ptype002'}; % auto-names used when N_i is numeric
V_final_PT=struct(); Policy_final_PT=struct();
for ii=1:N_i
    V_final_PT.(names_PT{ii})=V_final;
    Policy_final_PT.(names_PT{ii})=Policy_final;
end

%% PType solve
[VPath_PT,PolicyPath_PT]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, Policy_final_PT, Params, n_d, n_a, n_z, N_j, N_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn_PT, transpathoptionsbaseline, vfoptions);
Policy_init_PT=struct();
for ii=1:N_i
    Policy_init_PT.(names_PT{ii})=PolicyPath_PT.(names_PT{ii})(:,:,:,:,1);
end
AgentDist_initial_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_init_PT,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
AgentDistPath_PT=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist, PricePath, ParamPath, PolicyPath_PT, AgeWeightParamNames, PTypeDistParamNames,n_d,n_a,n_z,N_j,N_i,pi_z, T,Params, transpathoptionsbaseline, simoptions);
AggVarsPath_PT=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate_PT, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, N_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);

%% Two separate single-type solves
ptw=Params.(PTypeDistParamNames{:});
VPath_solo=cell(N_i,1); PolicyPath_solo=cell(N_i,1); AgentDistPath_solo=cell(N_i,1); AggVarsPath_solo=cell(N_i,1);
for ii=1:N_i
    Params_ii=Params;
    Params_ii.kappa_j=Params.kappa_j_pt(ii,:);
    [VPath_solo{ii},PolicyPath_solo{ii}]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params_ii, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn_NoPT, transpathoptionsbaseline, vfoptions);
    AgentDist_initial_ii=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyPath_solo{ii}(:,:,:,:,1),n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions);
    AgentDistPath_solo{ii}=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial_ii, jequaloneDist, PricePath, ParamPath, PolicyPath_solo{ii}, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params_ii, transpathoptionsbaseline, simoptions);
    AggVarsPath_solo{ii}=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate_NoPT, AgentDistPath_solo{ii}, PolicyPath_solo{ii}, PricePath, ParamPath, Params_ii, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
end

%% Compare per-type slices
for ii=1:N_i
    nm=names_PT{ii};
    fprintf('Two PTypes vs two solves, VPath        (type %d), this should be zero: %.3e \n',ii,max(abs(VPath_PT.(nm)(:)-VPath_solo{ii}(:))))
    fprintf('Two PTypes vs two solves, PolicyPath   (type %d), this should be zero: %.3e \n',ii,max(abs(PolicyPath_PT.(nm)(:)-PolicyPath_solo{ii}(:))))
    fprintf('Two PTypes vs two solves, AgentDistPath(type %d), this should be zero: %.3e \n',ii,max(abs(AgentDistPath_PT.(nm)(:)-AgentDistPath_solo{ii}(:))))
end

%% Aggregated AggVarsPath.Mean equals ptw-weighted combination of solo mean paths
agg_assets_mean  =ptw(1)*AggVarsPath_solo{1}.assets.Mean(:)  +ptw(2)*AggVarsPath_solo{2}.assets.Mean(:);
agg_earnings_mean=ptw(1)*AggVarsPath_solo{1}.earnings.Mean(:)+ptw(2)*AggVarsPath_solo{2}.earnings.Mean(:);
fprintf('Two PTypes vs two solves, AggVarsPath assets.Mean,   this should be zero: %.3e \n',max(abs(AggVarsPath_PT.assets.Mean(:)  -agg_assets_mean)))
fprintf('Two PTypes vs two solves, AggVarsPath earnings.Mean, this should be zero: %.3e \n',max(abs(AggVarsPath_PT.earnings.Mean(:)-agg_earnings_mean)))

%% Do some graphs of the AggVars path to see them
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:T,AggVarsPath_PT.earnings.Mean, 1:1:T,agg_earnings_mean)
title('Earnings Mean (PType vs weighted two solves)')
legend('PType','weighted solo')
subplot(2,1,2); plot(1:1:T,AggVarsPath_PT.assets.Mean, 1:1:T,agg_assets_mean)
title('Assets Mean (PType vs weighted two solves)')
legend('PType','weighted solo')

%% Statistics along the path, with a parameter that differs by ptype given as a matrix (kappa_j_pt, N_i-by-N_j)
%   (PA2) each ptype's AllStats and life-cycle profiles vs the commands without permanent types on its solo paths
%   (PA4) on a do-nothing path (prices and parameters held at their Params values) that starts from the
%         stationary PType agent dist, every period's AllStats is the stationary PType AllStats, pooled and by ptype
%   (PL3) on the same path, every cohort's profile is the stationary PType profile, where the cohort is alive
statsname='two types, kappa_j_pt';
OutputTypes={'Mean','Median','Variance','LorenzCurve','Gini','QuantileCutoffs','QuantileMeans'};
scalarstat=[1,1,1,0,1,0,0]; % Mean, Median, Variance, Gini are one number per age
fnnames=fieldnames(FnsToEvaluate_PT);
AllStats_PT=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsToEvaluate_PT, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, N_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
LCP_PT=LifeCycleProfiles_TransPath_FHorz_Case1_PType(FnsToEvaluate_PT, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, N_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
for ii=1:N_i
    nm=names_PT{ii};
    Params_ii=Params;
    Params_ii.kappa_j=Params.kappa_j_pt(ii,:);
    REF=EvalFnOnTransPath_AllStats_Case1_FHorz(FnsToEvaluate_NoPT, AgentDistPath_solo{ii}, PolicyPath_solo{ii}, PricePath, ParamPath, Params_ii, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
    REFL=LifeCycleProfiles_TransPath_FHorz_Case1(FnsToEvaluate_NoPT, AgentDistPath_solo{ii}, PolicyPath_solo{ii}, PricePath, ParamPath, Params_ii, T, n_d, n_a, n_z, N_j, d_grid, a_grid, z_grid, transpathoptionsbaseline, simoptions);
    % every field of the reference, to any depth, found in the PType output with the ptype name after the FnsToEvaluate name
    maxdiff=0; nanmismatch=0;
    queue={};
    REFnames=fieldnames(REF);
    for kk=1:length(REFnames)
        queue{end+1}=REFnames(kk);
    end
    while ~isempty(queue)
        fieldpath=queue{1};
        queue(1)=[];
        O=getfield(REF,fieldpath{:});
        if isstruct(O)
            Onames=fieldnames(O);
            for kk=1:length(Onames)
                queue{end+1}=[fieldpath,Onames(kk)];
            end
        else
            Ppath=[fieldpath(1),{nm},fieldpath(2:end)];
            P=AllStats_PT;
            found=1;
            for kk=1:length(Ppath)
                if isstruct(P) && isfield(P,Ppath{kk})
                    P=P.(Ppath{kk});
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
    fprintf('PType stats, %s: (PA2) %s vs its solo paths, EvalFnOnTransPath_AllStats_PType vs without permanent types, every stat, this should be zero: %.3e \n',statsname,nm,maxdiff)
    fprintf('PType stats, %s: (PA2) %s vs its solo paths, EvalFnOnTransPath_AllStats_PType vs without permanent types, every stat, positions NaN (or text) in one but not the other, or missing, this should be zero: %i \n',statsname,nm,nanmismatch)
    maxdiff=0; nanmismatch=0;
    for ff=1:length(fnnames)
        for oo=1:length(OutputTypes)
            PI=gather(LCP_PT.(fnnames{ff}).(nm).initialages.(OutputTypes{oo}));
            PB=gather(LCP_PT.(fnnames{ff}).(nm).bornduringtranstion.(OutputTypes{oo}));
            OI=gather(REFL.(fnnames{ff}).initialages.(OutputTypes{oo}));
            OB=gather(REFL.(fnnames{ff}).bornduringtranstion.(OutputTypes{oo}));
            maxdiff=max([maxdiff; abs(PI(:)-OI(:)); abs(PB(:)-OB(:))],[],'omitnan');
            nanmismatch=nanmismatch+sum(isnan(PI(:))~=isnan(OI(:)))+sum(isnan(PB(:))~=isnan(OB(:)));
        end
    end
    fprintf('PType stats, %s: (PA2) %s vs its solo paths, LifeCycleProfiles_TransPath_PType vs without permanent types, all seven outputs, this should be zero: %.3e \n',statsname,nm,maxdiff)
    fprintf('PType stats, %s: (PA2) %s vs its solo paths, LifeCycleProfiles_TransPath_PType vs without permanent types, entries NaN in one but not the other (the alive pattern), this should be zero: %i \n',statsname,nm,nanmismatch)
end

% (PA4)/(PL3) the do-nothing path
[V_s,Policy_s]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions);
AgentDist_s=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_s,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
PricePathNamesS=fieldnames(PricePath);
ParamPathNamesS=fieldnames(ParamPath);
PricePathC=struct();
for kk=1:length(PricePathNamesS)
    PricePathC.(PricePathNamesS{kk})=Params.(PricePathNamesS{kk})*ones(1,T);
end
ParamPathC=struct();
for kk=1:length(ParamPathNamesS)
    ParamPathC.(ParamPathNamesS{kk})=Params.(ParamPathNamesS{kk})*ones(1,T);
end
[~,PolicyPathC]=ValueFnOnTransPath_Case1_FHorz_PType(PricePathC, ParamPathC, T, V_s, Policy_s, Params, n_d, n_a, n_z, N_j, N_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn_PT, transpathoptionsbaseline, vfoptions);
AgentDistPathC=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_s, jequaloneDist, PricePathC, ParamPathC, PolicyPathC, AgeWeightParamNames, PTypeDistParamNames,n_d,n_a,n_z,N_j,N_i,pi_z, T,Params, transpathoptionsbaseline, simoptions);
AllStatsPathC=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsToEvaluate_PT, AgentDistPathC, PolicyPathC, PricePathC, ParamPathC, Params, T, n_d, n_a, n_z, N_j, N_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
LCPPathC=LifeCycleProfiles_TransPath_FHorz_Case1_PType(FnsToEvaluate_PT, AgentDistPathC, PolicyPathC, PricePathC, ParamPathC, Params, T, n_d, n_a, n_z, N_j, N_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
AS_stat=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(AgentDist_s,Policy_s,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);
LCP_stat=LifeCycleProfiles_FHorz_Case1_PType(AgentDist_s,Policy_s,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);
maxdiff=0; nanmismatch=0;
for tt=1:T
    queue={{}};
    while ~isempty(queue)
        fieldpath=queue{1};
        queue(1)=[];
        if isempty(fieldpath)
            O=AS_stat;
        else
            O=getfield(AS_stat,fieldpath{:});
        end
        if isstruct(O)
            Onames=fieldnames(O);
            for kk=1:length(Onames)
                queue{end+1}=[fieldpath,Onames(kk)];
            end
        elseif ischar(O) || isstring(O) || iscell(O) % text (the comments): must be the same text
            P=getfield(AllStatsPathC,fieldpath{:});
            nanmismatch=nanmismatch+(~isequal(P{tt},O));
        else
            O=gather(O);
            P=getfield(AllStatsPathC,fieldpath{:});
            P=P(:,tt);
            maxdiff=max([maxdiff; abs(P(:)-O(:))],[],'omitnan');
            nanmismatch=nanmismatch+sum(isnan(P(:))~=isnan(O(:)));
        end
    end
end
fprintf('PType stats, %s: (PA4) on the do-nothing path, EvalFnOnTransPath_AllStats_PType vs the stationary PType AllStats, every period, pooled and by ptype, this should be zero: %.3e \n',statsname,maxdiff)
fprintf('PType stats, %s: (PA4) on the do-nothing path, EvalFnOnTransPath_AllStats_PType vs the stationary PType AllStats, every period, positions NaN (or text) in one but not the other, this should be zero: %i \n',statsname,nanmismatch)
groupnames=[{''},names_PT];
maxdiff=0; nanmismatch=0;
for gg=1:length(groupnames)
    for ff=1:length(fnnames)
        if isempty(groupnames{gg})
            OUTff=LCPPathC.(fnnames{ff});
            temp_ff=LCP_stat.(fnnames{ff});
        else
            OUTff=LCPPathC.(fnnames{ff}).(groupnames{gg});
            temp_ff=LCP_stat.(fnnames{ff}).(groupnames{gg});
        end
        for oo=1:length(OutputTypes)
            PI=gather(OUTff.initialages.(OutputTypes{oo}));
            PB=gather(OUTff.bornduringtranstion.(OutputTypes{oo}));
            OI=nan(size(PI));
            OB=nan(size(PB));
            temp=gather(temp_ff.(OutputTypes{oo}));
            for tt=1:min(T,N_j)
                for jj=1:N_j-tt+1
                    if scalarstat(oo)==1
                        OI(jj,jj+tt-1)=temp(jj+tt-1);
                    else
                        OI(:,jj,jj+tt-1)=temp(:,jj+tt-1);
                    end
                end
            end
            for bb=1:T
                for age=1:N_j
                    if bb+age-1<=T
                        if scalarstat(oo)==1
                            OB(bb,age)=temp(age);
                        else
                            OB(:,bb,age)=temp(:,age);
                        end
                    end
                end
            end
            maxdiff=max([maxdiff; abs(PI(:)-OI(:)); abs(PB(:)-OB(:))],[],'omitnan');
            nanmismatch=nanmismatch+sum(isnan(PI(:))~=isnan(OI(:)))+sum(isnan(PB(:))~=isnan(OB(:)));
        end
    end
end
fprintf('PType stats, %s: (PL3) on the do-nothing path, LifeCycleProfiles_TransPath_PType vs the stationary PType profiles, pooled and by ptype, all seven outputs, this should be zero: %.3e \n',statsname,maxdiff)
fprintf('PType stats, %s: (PL3) on the do-nothing path, LifeCycleProfiles_TransPath_PType vs the stationary PType profiles, entries NaN in one but not the other (the alive pattern), this should be zero: %i \n',statsname,nanmismatch)

output=struct();

end
