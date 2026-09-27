function output=CoreFHorzPType_ValuesOnGridPerTypeShocks(n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames)
% Test: EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType when n_z differs by ptype (n_z a struct):
% one type with z, one without. Per-type ValuesOnGrid must equal solo ValuesOnGrid.
%
% Until 2026-09-27 the command did N_a=prod(n_a); N_z=prod(n_z) at the top, before the per-type
% loop, and prod of a struct is an error. Those two lines are only needed to preallocate for a
% non-struct FnsToEvaluate, and were moved inside that branch. This part is kept last so that,
% were it to error again, it could not stop any other part from running.

Names_i={'z','noz'};
N_i=length(Names_i);

n_d=0;
d_grid=[];

n_z_PT.z=n_z;      z_grid_PT.z=z_grid; pi_z_PT.z=pi_z;
n_z_PT.noz=0;      z_grid_PT.noz=[];   pi_z_PT.noz=[];

ReturnFn_PT.z  =@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
ReturnFn_PT.noz=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate_PT.assets.z  =@(aprime,a,z) a;
FnsToEvaluate_PT.assets.noz=@(aprime,a) a;
FnsToEvaluate_PT.income.z  =@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate_PT.income.noz=@(aprime,a,w,kappa_j) w*kappa_j;
FnNames=fieldnames(FnsToEvaluate_PT);

vfoptions=struct();
simoptions=struct();

[~,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z_PT,N_j,Names_i,d_grid,a_grid,z_grid_PT,pi_z_PT,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions);
ValuesOnGrid_PT=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z_PT,N_j,Names_i,d_grid,a_grid,z_grid_PT,simoptions);

for ii=1:N_i
    nm=Names_i{ii};
    [~,Policy_ii]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z_PT.(nm),N_j,d_grid,a_grid,z_grid_PT.(nm),pi_z_PT.(nm),ReturnFn_PT.(nm),Params,DiscountFactorParamNames,[],vfoptions);
    FnsToEvaluate_solo=struct();
    for ff=1:length(FnNames)
        FnsToEvaluate_solo.(FnNames{ff})=FnsToEvaluate_PT.(FnNames{ff}).(nm);
    end
    ValuesOnGrid_ii=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy_ii,FnsToEvaluate_solo,Params,[],n_d,n_a,n_z_PT.(nm),N_j,d_grid,a_grid,z_grid_PT.(nm),simoptions);
    for ff=1:length(FnNames)
        fn=FnNames{ff};
        fprintf('ValuesOnGrid per-type n_z, %s (type %s), this should be zero: %.3e \n',fn,nm,max(abs(ValuesOnGrid_PT.(fn).(nm)(:)-ValuesOnGrid_ii.(fn)(:))))
    end
end

output=struct();

end
