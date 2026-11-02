function output=CoreFHorzPType_NoApplicableFns(n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames)
% Test: a ptype that NONE of the FnsToEvaluate apply to. Until 2026-11-01 the PType commands crashed on it:
% PType_FnsToEvaluate returns an empty cell for such a ptype, and the callers passed it on (fieldnames({}), or a
% deterministic command that never assigns its function list when there are no functions). Found by AACP2022, whose
% FnsToEvaluate_Cpath has household entries only. The callers now skip such a ptype (or, in the general eqm commands,
% skip its aggregate evaluation).
% THE ORACLE: ABSENT EQUALS ZERO. The same FnsToEvaluate with the function DEFINED for the other type as zero never
% takes the skip, and every pooled output must be the same: a type with no applicable function adds nothing.
%   EvalFnOnAgentDist_AggVars_FHorz_Case1_PType      pooled Mean identical; the type with no function has Mean zero
%   EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType the low type's values identical; no field for the high type
%   SimPanelValues_FHorz_Case1_PType                 the high type's simulations NaN; and a MIXED case (one function of
%       both types, one of the low type only), which used to crash on the row index (WhichFnsForCurrentPType is 0 for
%       a function that does not apply, and was used to index the rows) even though no type had zero functions
%   HeteroAgentStationaryEqm_Case1_FHorz_PType       the general eqm (an economy where only the low type supplies
%       capital and labour) identical with the high type's functions absent and zero, in p_eqm and the conditions; and
%       with a condition conditional on ptype (heteroagentoptions.GEptype, the GEptype subfn) at the initial prices
% The two types are those of part 4 (per-type kappa_j_pt), no d, with z.

Names_i={'low','high'};
n_d=0; d_grid=[];

ReturnFn=@(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension);

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

vfoptions=struct();
simoptions=struct();

[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);

% Absent: the high type has no function at all. Zero: the same, with the high type's functions defined as zero.
FnsAbsent=struct();
FnsAbsent.assetslow.low=@(aprime,a,z) a;
FnsAbsent.earnlow.low=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsZero=FnsAbsent;
FnsZero.assetslow.high=@(aprime,a,z) 0*a;
FnsZero.earnlow.high=@(aprime,a,z,w,kappa_j_pt) 0*z;
FnNames=fieldnames(FnsAbsent);

%% AggVars
AggAbsent=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist,Policy,FnsAbsent,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
AggZero=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist,Policy,FnsZero,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
for ff=1:length(FnNames)
    fn=FnNames{ff};
    fprintf('No applicable fns, AggVars %s: pooled Mean with the high type absent vs zero, this should be zero: %.3e \n',fn,abs(AggAbsent.(fn).Mean-AggZero.(fn).Mean))
    fprintf('No applicable fns, AggVars %s: and the low type''s own Mean, this should be zero: %.3e \n',fn,abs(AggAbsent.(fn).low.Mean-AggZero.(fn).low.Mean))
    fprintf('No applicable fns, AggVars %s: the high type''s Mean (zero, as its column of AggVarsFull is never filled), this should be zero: %.3e \n',fn,abs(AggAbsent.(fn).high.Mean))
end
fprintf('No applicable fns, AggVars: for scale, the low type''s assets are not zero, should NOT be zero: %.3e \n',abs(AggAbsent.assetslow.Mean))

%% ValuesOnGrid
VoGAbsent=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(Policy,FnsAbsent,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
VoGZero=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(Policy,FnsZero,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
for ff=1:length(FnNames)
    fn=FnNames{ff};
    fprintf('No applicable fns, ValuesOnGrid %s: the low type with the high type absent vs zero, this should be zero: %.3e \n',fn,max(abs(gather(VoGAbsent.(fn).low(:))-gather(VoGZero.(fn).low(:)))))
    fprintf('No applicable fns, ValuesOnGrid %s: the high type has no field, should be one: %i \n',fn,~isfield(VoGAbsent.(fn),'high'))
end

%% SimPanelValues: the high type absent, and the mixed case
% PType_numbersims=floor(ptypeweights*numbersims): the first ns(1) simulations are the low type, the next ns(2) the high type
simoptions_panel=simoptions;
simoptions_panel.numbersims=1000;
ns=floor(Params.(PTypeDistParamNames{1})*simoptions_panel.numbersims);
lowcols=1:ns(1);
highcols=ns(1)+1:ns(1)+ns(2);
PanelAbsent=SimPanelValues_FHorz_Case1_PType(jequaloneDist,PTypeDistParamNames,Policy,FnsAbsent,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_panel);
fprintf('No applicable fns, SimPanel with the high type absent: the low type''s simulations are finite, should be one: %i \n',all(isfinite(PanelAbsent.assetslow(:,lowcols)),'all'))
fprintf('No applicable fns, SimPanel with the high type absent: the high type''s simulations are NaN, should be one: %i \n',all(isnan(PanelAbsent.assetslow(:,highcols)),'all'))
fprintf('No applicable fns, SimPanel with the high type absent: and of the second function too (the panel was set up with one row, length of a struct, and the rows it grew were zero, not NaN), should be one: %i \n',all(isnan(PanelAbsent.earnlow(:,highcols)),'all'))
FnsMixed=struct();
FnsMixed.assets=@(aprime,a,z) a; % both types
FnsMixed.earnlow.low=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z; % the low type only
PanelMixed=SimPanelValues_FHorz_Case1_PType(jequaloneDist,PTypeDistParamNames,Policy,FnsMixed,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_panel);
fprintf('No applicable fns, SimPanel mixed: the function of both types is finite for both types, should be one: %i \n',all(isfinite(PanelMixed.assets(:,[lowcols,highcols])),'all'))
fprintf('No applicable fns, SimPanel mixed: the low-type function is finite for the low type, should be one: %i \n',all(isfinite(PanelMixed.earnlow(:,lowcols)),'all'))
fprintf('No applicable fns, SimPanel mixed: and NaN for the high type, should be one: %i \n',all(isnan(PanelMixed.earnlow(:,highcols)),'all'))

%% AllStats: the grouped Mean and StdDeviation are over the types a function applies to
% With the FIRST type absent (the function only of the high type) the grouped StdDeviation's cross term uses the
% difference of the Means of the earlier types, which was NaN for an absent type (as was the grouped Mean, with
% either type absent, through 0*NaN). Only one type contributes, so the grouped stats equal its own.
FnsHighOnly=struct();
FnsHighOnly.assetshigh.high=@(aprime,a,z) a;
AllStatsLowOnly=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsAbsent,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
AllStatsHighOnly=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsHighOnly,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
fprintf('No applicable fns, AllStats, the second type absent: grouped Mean vs the low type''s, this should be zero: %.3e \n',abs(AllStatsLowOnly.assetslow.Mean-AllStatsLowOnly.assetslow.low.Mean))
fprintf('No applicable fns, AllStats, the second type absent: grouped StdDeviation vs the low type''s, this should be zero: %.3e \n',abs(AllStatsLowOnly.assetslow.StdDeviation-AllStatsLowOnly.assetslow.low.StdDeviation))
fprintf('No applicable fns, AllStats, the first type absent: grouped Mean vs the high type''s, this should be zero: %.3e \n',abs(AllStatsHighOnly.assetshigh.Mean-AllStatsHighOnly.assetshigh.high.Mean))
fprintf('No applicable fns, AllStats, the first type absent: grouped StdDeviation vs the high type''s, this should be zero: %.3e \n',abs(AllStatsHighOnly.assetshigh.StdDeviation-AllStatsHighOnly.assetshigh.high.StdDeviation))
fprintf('No applicable fns, AllStats: for scale, the high type''s StdDeviation of assets is not zero, should NOT be zero: %.3e \n',AllStatsHighOnly.assetshigh.high.StdDeviation)

%% General eqm: only the low type supplies capital and labour
ParamsGE=Params;
ParamsGE.firmalpha=0.36;
ParamsGE.firmdelta=0.05;
GeneralEqmEqns.CapitalMarket=@(r,K,L,firmalpha,firmdelta) r-(firmalpha*(K^(firmalpha-1))*(L^(1-firmalpha))-firmdelta);
GeneralEqmEqns.LabourMarket=@(w,K,L,firmalpha) w-((1-firmalpha)*(K^firmalpha)*(L^(-firmalpha)));
FnsGEAbsent=struct();
FnsGEAbsent.K.low=@(aprime,a,z) a;
FnsGEAbsent.L.low=@(aprime,a,z,kappa_j_pt) kappa_j_pt*z;
FnsGEZero=FnsGEAbsent;
FnsGEZero.K.high=@(aprime,a,z) 0*a;
FnsGEZero.L.high=@(aprime,a,z,kappa_j_pt) 0*z;
% The oracle is an identity between the two solves, which follow the same arithmetic whether or not the solver
% converges, so a short solve is enough
heteroagentoptions=struct();
heteroagentoptions.maxiter=30;
[p_Absent,GEc_Absent]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,0,pi_z,d_grid,a_grid,z_grid,jequaloneDist,ReturnFn,FnsGEAbsent,GeneralEqmEqns,ParamsGE,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,{'r','w'},heteroagentoptions,simoptions,vfoptions);
[p_Zero,GEc_Zero]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,0,pi_z,d_grid,a_grid,z_grid,jequaloneDist,ReturnFn,FnsGEZero,GeneralEqmEqns,ParamsGE,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,{'r','w'},heteroagentoptions,simoptions,vfoptions);
fprintf('No applicable fns, general eqm: p_eqm with the high type''s functions absent vs zero, this should be zero: %.3e \n',max(abs([p_Absent.r-p_Zero.r,p_Absent.w-p_Zero.w])))
fprintf('No applicable fns, general eqm: and the general eqm conditions, this should be zero: %.3e \n',max(abs([GEc_Absent.CapitalMarket-GEc_Zero.CapitalMarket,GEc_Absent.LabourMarket-GEc_Zero.LabourMarket])))
fprintf('No applicable fns, general eqm: for information, r=%.6f, w=%.6f, conditions %.3e %.3e \n',p_Absent.r,p_Absent.w,GEc_Absent.CapitalMarket,GEc_Absent.LabourMarket)

%% A general eqm condition conditional on ptype (HeteroAgentStationaryEqm_Case1_FHorz_PType_GEptype_subfn), evaluated at
% the initial prices (heteroagentoptions.maxiter=0). KReport is by ptype: each type's own K against Kbar.
ParamsGE.Kbar=1;
GeneralEqmEqnsPT=GeneralEqmEqns;
GeneralEqmEqnsPT.KReport=@(K,Kbar) K-Kbar;
heteroagentoptionsPT=struct();
heteroagentoptionsPT.maxiter=0;
heteroagentoptionsPT.GEptype={'KReport'};
[~,GEcPT_Absent]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,0,pi_z,d_grid,a_grid,z_grid,jequaloneDist,ReturnFn,FnsGEAbsent,GeneralEqmEqnsPT,ParamsGE,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,{'r','w'},heteroagentoptionsPT,simoptions,vfoptions);
[~,GEcPT_Zero]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,0,pi_z,d_grid,a_grid,z_grid,jequaloneDist,ReturnFn,FnsGEZero,GeneralEqmEqnsPT,ParamsGE,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,{'r','w'},heteroagentoptionsPT,simoptions,vfoptions);
fprintf('No applicable fns, GEptype: the economy-wide conditions with the high type''s functions absent vs zero, this should be zero: %.3e \n',max(abs([GEcPT_Absent.CapitalMarket-GEcPT_Zero.CapitalMarket,GEcPT_Absent.LabourMarket-GEcPT_Zero.LabourMarket])))
fprintf('No applicable fns, GEptype: the by-ptype condition (both types), this should be zero: %.3e \n',max(abs(GEcPT_Absent.KReport(:)-GEcPT_Zero.KReport(:))))
fprintf('No applicable fns, GEptype: the by-ptype condition has one value per type, should be one: %i \n',numel(GEcPT_Absent.KReport)==2)
fprintf('No applicable fns, GEptype: the high type''s K is zero, so its condition is -Kbar, this should be zero: %.3e \n',abs(GEcPT_Absent.KReport(2)+ParamsGE.Kbar))
fprintf('No applicable fns, GEptype: the low type''s K is not zero, should NOT be zero: %.3e \n',abs(GEcPT_Absent.KReport(1)+ParamsGE.Kbar))

output=struct();

end
