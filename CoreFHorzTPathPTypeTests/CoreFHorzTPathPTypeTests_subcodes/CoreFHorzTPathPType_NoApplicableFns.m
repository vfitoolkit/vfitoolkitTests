function output=CoreFHorzTPathPType_NoApplicableFns(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,transpathoptionsbaseline)
% Test: a ptype that NONE of the FnsToEvaluate apply to, in the FHorz PType commands on the transition path. Until
% 2026-11-01 they crashed on it (PType_FnsToEvaluate returns an empty cell for such a ptype). Found by AACP2022, whose
% FnsToEvaluate_Cpath has household entries only. Now EvalFnOnTransPath_AggVars skips such a ptype, and AllStats and
% LifeCycleProfiles_TransPath give it an empty FnsToEvaluate for its lags (it still needs its lagsinfo).
% Two genuinely different types (per-type kappa_j_pt), no d, with z. FnsAbsent apply to the low type only; FnsZero are
% the same plus the high type's functions defined as zero.
%   EvalFnOnTransPath_AggVars_Case1_FHorz_PType: ABSENT EQUALS ZERO. The pooled Mean is the ptypeweights-weighted sum of
%       the per-type Means, so a zero function adds nothing: pooled and low Means identical; no high field when absent.
%   EvalFnOnTransPath_AllStats_Case1_FHorz_PType and LifeCycleProfiles_TransPath_FHorz_Case1_PType: the pooled stats are
%       over the types the function applies to (the masses renormalised by SigmaNxi), so absent is NOT zero for them
%       (a zero type is a mass at zero). The oracle is instead: every low-type stat identical with the high type absent
%       vs zero; no high field (or only NaN) when absent; and the pooled Mean and StdDeviation (Variance for LifeCycleProfiles_TransPath) with the high type absent
%       equal the low type's own (only the low type contributes).

Names_i={'low','high'};
N_i=2;
n_d=0;
d_grid=[];

ReturnFn=@(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension);

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

vfoptions=struct();
simoptions=struct();

% Trivial period-T V and Policy (as in CoreFHorzTPathPType_PerTypeFnsToEvaluate)
V_final_PT=struct(); Policy_final_PT=struct();
for ii=1:N_i
    V_final_PT.(Names_i{ii})=zeros([n_a,n_z,N_j],'gpuArray');
    Policy_final_PT.(Names_i{ii})=ones([1,n_a,n_z,N_j],'gpuArray');
end

[~,PolicyPath_PT]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, Policy_final_PT, Params, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
Policy_init_PT=struct();
for ii=1:N_i
    Policy_init_PT.(Names_i{ii})=PolicyPath_PT.(Names_i{ii})(:,:,:,:,1);
end
AgentDist_initial_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_init_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AgentDistPath_PT=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist, PricePath, ParamPath, PolicyPath_PT, AgeWeightParamNames, PTypeDistParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z, T,Params, transpathoptionsbaseline, simoptions);

% Absent: the high type has no function at all. Zero: the same, with the high type's functions defined as zero.
FnsAbsent=struct();
FnsAbsent.assetslow.low=@(aprime,a,z) a;
FnsAbsent.earnlow.low=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsZero=FnsAbsent;
FnsZero.assetslow.high=@(aprime,a,z) 0*a;
FnsZero.earnlow.high=@(aprime,a,z,w,kappa_j_pt) 0*z;
FnNames=fieldnames(FnsAbsent);

%% AggVars on the transition path: absent equals zero
AggAbsent=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsAbsent, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
AggZero=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsZero, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
for ff=1:length(FnNames)
    fn=FnNames{ff};
    fprintf('No applicable fns, AggVarsPath %s: pooled Mean with the high type absent vs zero, this should be zero: %.3e \n',fn,max(abs(AggAbsent.(fn).Mean(:)-AggZero.(fn).Mean(:))))
    fprintf('No applicable fns, AggVarsPath %s: and the low type''s own Mean, this should be zero: %.3e \n',fn,max(abs(AggAbsent.(fn).low.Mean(:)-AggZero.(fn).low.Mean(:))))
    fprintf('No applicable fns, AggVarsPath %s: the high type has no field, should be one: %i \n',fn,~isfield(AggAbsent.(fn),'high'))
end
fprintf('No applicable fns, AggVarsPath: for scale, the low type''s assets are not zero, should NOT be zero: %.3e \n',max(abs(AggAbsent.assetslow.Mean(:))))

%% AllStats and life-cycle profiles on the transition path
AllStatsAbsent=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsAbsent, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
AllStatsZero=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsZero, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
LCPAbsent=LifeCycleProfiles_TransPath_FHorz_Case1_PType(FnsAbsent, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);
LCPZero=LifeCycleProfiles_TransPath_FHorz_Case1_PType(FnsZero, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptionsbaseline, simoptions);

% Every low-type stat (to any depth) with the high type absent must be in the zero output with the same value
cmpcases={'AllStatsPath', AllStatsAbsent, AllStatsZero;
          'AgeCondStatsPath', LCPAbsent, LCPZero};
for cc=1:size(cmpcases,1)
    for ff=1:length(FnNames)
        fn=FnNames{ff};
        REF=cmpcases{cc,2}.(fn).low;
        OUT=cmpcases{cc,3}.(fn).low;
        maxdiff=0; nanmismatch=0; nleaves=0;
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
                nleaves=nleaves+1;
                if found==0
                    nanmismatch=nanmismatch+1;
                elseif iscell(O) || ischar(O) || isstring(O) % text (the comments): must be the same text
                    nanmismatch=nanmismatch+(~isequal(P,O));
                else
                    P=gather(P); O=gather(O);
                    maxdiff=max([maxdiff; abs(P(:)-O(:))],[],'omitnan');
                    nanmismatch=nanmismatch+sum(isnan(P(:))~=isnan(O(:)));
                end
            end
        end
        fprintf('No applicable fns, %s %s: every low-type stat with the high type absent vs zero (%i leaves), this should be zero: %.3e \n',cmpcases{cc,1},fn,nleaves,maxdiff)
        fprintf('No applicable fns, %s %s: low-type positions NaN (or text) in one but not the other, or missing, this should be zero: %i \n',cmpcases{cc,1},fn,nanmismatch)
        A=cmpcases{cc,2}.(fn);
        if isfield(A,'high')
            if strcmp(cmpcases{cc,1},'AllStatsPath')
                highvals=A.high.Mean;
            else
                highvals=[A.high.initialages.Mean(:); A.high.bornduringtranstion.Mean(:)];
            end
            nohigh=all(isnan(gather(highvals(:))));
        else
            nohigh=1;
        end
        fprintf('No applicable fns, %s %s: the high type has no field (or only NaN), should be one: %i \n',cmpcases{cc,1},fn,nohigh)
    end
end

% The pooled stats with the high type absent are over the low type alone, so its Mean and StdDeviation equal the low type's
for ff=1:length(FnNames)
    fn=FnNames{ff};
    A=AllStatsAbsent.(fn);
    fprintf('No applicable fns, AllStatsPath %s: pooled Mean with the high type absent vs the low type''s, this should be zero: %.3e \n',fn,max(abs(gather(A.Mean(:))-gather(A.low.Mean(:)))))
    fprintf('No applicable fns, AllStatsPath %s: pooled StdDeviation with the high type absent vs the low type''s, this should be zero: %.3e \n',fn,max(abs(gather(A.StdDeviation(:))-gather(A.low.StdDeviation(:)))))
    A=LCPAbsent.(fn);
    PA=gather([A.initialages.Mean(:); A.bornduringtranstion.Mean(:)]);
    PL=gather([A.low.initialages.Mean(:); A.low.bornduringtranstion.Mean(:)]);
    fprintf('No applicable fns, AgeCondStatsPath %s: pooled Mean with the high type absent vs the low type''s, this should be zero: %.3e \n',fn,max(abs(PA-PL),[],'omitnan'))
    fprintf('No applicable fns, AgeCondStatsPath %s: pooled Mean NaN where the low type''s is not, or vice versa, this should be zero: %i \n',fn,sum(isnan(PA)~=isnan(PL)))
    PA=gather([A.initialages.Variance(:); A.bornduringtranstion.Variance(:)]);
    PL=gather([A.low.initialages.Variance(:); A.low.bornduringtranstion.Variance(:)]);
    fprintf('No applicable fns, AgeCondStatsPath %s: pooled Variance (LifeCycleProfiles_TransPath reports Variance, not StdDeviation) with the high type absent vs the low type''s, this should be zero: %.3e \n',fn,max(abs(PA-PL),[],'omitnan'))
end
fprintf('No applicable fns, AllStatsPath: for scale, the low type''s mean assets are not zero, should NOT be zero: %.3e \n',max(abs(gather(AllStatsAbsent.assetslow.Mean(:)))))

output=struct();

end
