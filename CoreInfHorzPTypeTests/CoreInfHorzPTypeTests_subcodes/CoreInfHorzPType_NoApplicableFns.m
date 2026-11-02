function output=CoreInfHorzPType_NoApplicableFns(n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT)
% Test: the InfHorz PType general eqm with a ptype that NONE of the FnsToEvaluate apply to. Until 2026-11-01 the general
% eqm subfns passed such a ptype's empty FnsToEvaluate on to EvalFnOnAgentDist_AggVars_InfHorz, which crashed (found by
% AACP2022, whose FnsToEvaluate_Cpath has household entries only). The subfns now skip that ptype's aggregate evaluation.
% THE ORACLE: ABSENT EQUALS ZERO. FnsAbsent apply to the low type only; FnsZero are the same plus the high type's
% functions defined as zero. The aggregates are ptypeweights-weighted sums of the per-type ones (and, with GEptype, the
% per-type ones themselves), so a zero function adds nothing and the two must agree exactly.
%   HeteroAgentStationaryEqm_InfHorz_PType_subfn         a short solve of an economy where only the low type supplies
%       capital: p_eqm and the general eqm conditions identical
%   HeteroAgentStationaryEqm_InfHorz_PType_GEptype_subfn the general eqm conditions at the initial prices
%       (heteroagentoptions.maxiter=0) with a condition conditional on ptype (heteroagentoptions.GEptype): identical,
%       and the high type's by-ptype condition sees zero capital
%   EvalFnOnAgentDist_AllStats_InfHorz_PType            the grouped Mean and StdDeviation with one type absent (either
%       order) equal the other type's own
% Two types differing by the per-type wage w_pt of the bank setup, no d, with z.

Names_i={'low','high'};
n_d=0; d_grid=[];
n_p=0;

ReturnFn=@(aprime,a,z,r,w_pt,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w_pt,sigma);

ParamsGE=Params;
ParamsGE.firmalpha=0.36;
ParamsGE.firmdelta=0.05;
ParamsGE.Lbar=1; % labour is exogenous, so r is the only general eqm price
ParamsGE.Kbar=1; % a reference level for the by-ptype condition

FnsAbsent=struct();
FnsAbsent.K.low=@(aprime,a,z) a;
FnsZero=FnsAbsent;
FnsZero.K.high=@(aprime,a,z) 0*a;

GeneralEqmEqns=struct();
GeneralEqmEqns.CapitalMarket=@(r,K,Lbar,firmalpha,firmdelta) r-(firmalpha*(K^(firmalpha-1))*(Lbar^(1-firmalpha))-firmdelta);

vfoptions=struct();

%% The pooled general eqm (HeteroAgentStationaryEqm_InfHorz_PType_subfn): a short solve
% The oracle is an identity between the two solves, which follow the same arithmetic whether or not the solver
% converges, so a short solve is enough
heteroagentoptions=struct();
heteroagentoptions.maxiter=30;
[p_Absent,GEc_Absent]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsAbsent, GeneralEqmEqns, ParamsGE, DiscountFactorParamNames, PTypeDistParamNames, {'r'},heteroagentoptions, simoptionsPT, vfoptions);
[p_Zero,GEc_Zero]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsZero, GeneralEqmEqns, ParamsGE, DiscountFactorParamNames, PTypeDistParamNames, {'r'},heteroagentoptions, simoptionsPT, vfoptions);
fprintf('No applicable fns, InfHorz general eqm: p_eqm with the high type''s function absent vs zero, this should be zero: %.3e \n',abs(p_Absent.r-p_Zero.r))
fprintf('No applicable fns, InfHorz general eqm: and the general eqm condition, this should be zero: %.3e \n',abs(GEc_Absent.CapitalMarket-GEc_Zero.CapitalMarket))
fprintf('No applicable fns, InfHorz general eqm: for scale, the solve moved r from %.6f, should NOT be zero: %.3e \n',ParamsGE.r,abs(p_Absent.r-ParamsGE.r))
fprintf('No applicable fns, InfHorz general eqm: for information, r=%.6f, condition %.3e \n',p_Absent.r,GEc_Absent.CapitalMarket)

%% A general eqm condition conditional on ptype (HeteroAgentStationaryEqm_InfHorz_PType_GEptype_subfn), evaluated
% at the initial prices. KReport is by ptype: each type's own K against Kbar.
GeneralEqmEqnsPT=GeneralEqmEqns;
GeneralEqmEqnsPT.KReport=@(K,Kbar) K-Kbar;
heteroagentoptionsPT=struct();
heteroagentoptionsPT.maxiter=0;
heteroagentoptionsPT.GEptype={'KReport'};
[~,GEcPT_Absent]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsAbsent, GeneralEqmEqnsPT, ParamsGE, DiscountFactorParamNames, PTypeDistParamNames, {'r'},heteroagentoptionsPT, simoptionsPT, vfoptions);
[~,GEcPT_Zero]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsZero, GeneralEqmEqnsPT, ParamsGE, DiscountFactorParamNames, PTypeDistParamNames, {'r'},heteroagentoptionsPT, simoptionsPT, vfoptions);
fprintf('No applicable fns, InfHorz GEptype: the economy-wide condition with the high type''s function absent vs zero, this should be zero: %.3e \n',abs(GEcPT_Absent.CapitalMarket-GEcPT_Zero.CapitalMarket))
fprintf('No applicable fns, InfHorz GEptype: the by-ptype condition (both types), this should be zero: %.3e \n',max(abs(GEcPT_Absent.KReport(:)-GEcPT_Zero.KReport(:))))
fprintf('No applicable fns, InfHorz GEptype: the by-ptype condition has one value per type, should be one: %i \n',numel(GEcPT_Absent.KReport)==2)
fprintf('No applicable fns, InfHorz GEptype: the high type''s K is zero, so its condition is -Kbar, this should be zero: %.3e \n',abs(GEcPT_Absent.KReport(2)+ParamsGE.Kbar))
fprintf('No applicable fns, InfHorz GEptype: the low type''s K is not zero, should NOT be zero: %.3e \n',abs(GEcPT_Absent.KReport(1)+ParamsGE.Kbar))

%% EvalFnOnAgentDist_AllStats_InfHorz_PType: the grouped Mean and StdDeviation are over the types a function applies to
% A ptype the function does not apply to has a NaN Mean, and the grouped Mean (sum of the indicator times the ptype
% weight times the Mean) was NaN through 0*NaN, with either type absent; the grouped StdDeviation too, through the
% difference of the Means of the earlier types, when an earlier type is absent (found 2026-11-03 on the FHorz path).
% Only one type contributes, so the grouped stats equal its own. Both orders: the second type absent, the first absent.
[~,Policy_PT]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
Dist_PT=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,Names_i,pi_z,Params,simoptionsPT);
FnsLowOnly=struct();
FnsLowOnly.assetslow.low=@(aprime,a,z) a;
FnsHighOnly=struct();
FnsHighOnly.assetshigh.high=@(aprime,a,z) a;
AllStatsLowOnly=EvalFnOnAgentDist_AllStats_InfHorz_PType(Dist_PT,Policy_PT,FnsLowOnly,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,struct());
AllStatsHighOnly=EvalFnOnAgentDist_AllStats_InfHorz_PType(Dist_PT,Policy_PT,FnsHighOnly,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,struct());
fprintf('No applicable fns, InfHorz AllStats, the second type absent: grouped Mean vs the low type''s, this should be zero: %.3e \n',abs(AllStatsLowOnly.assetslow.Mean-AllStatsLowOnly.assetslow.low.Mean))
fprintf('No applicable fns, InfHorz AllStats, the second type absent: grouped StdDeviation vs the low type''s, this should be zero: %.3e \n',abs(AllStatsLowOnly.assetslow.StdDeviation-AllStatsLowOnly.assetslow.low.StdDeviation))
fprintf('No applicable fns, InfHorz AllStats, the first type absent: grouped Mean vs the high type''s, this should be zero: %.3e \n',abs(AllStatsHighOnly.assetshigh.Mean-AllStatsHighOnly.assetshigh.high.Mean))
fprintf('No applicable fns, InfHorz AllStats, the first type absent: grouped StdDeviation vs the high type''s, this should be zero: %.3e \n',abs(AllStatsHighOnly.assetshigh.StdDeviation-AllStatsHighOnly.assetshigh.high.StdDeviation))
fprintf('No applicable fns, InfHorz AllStats: for scale, the high type''s StdDeviation of assets is not zero, should NOT be zero: %.3e \n',AllStatsHighOnly.assetshigh.high.StdDeviation)

output=struct();

end
