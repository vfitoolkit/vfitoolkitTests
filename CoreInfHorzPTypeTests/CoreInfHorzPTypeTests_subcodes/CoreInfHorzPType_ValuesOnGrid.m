function output=CoreInfHorzPType_ValuesOnGrid(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT)
% Test: EvalFnOnAgentDist_ValuesOnGrid_InfHorz_PType (new): ValuesOnGrid.(fnname).(typename) is the single-type
% EvalFnOnAgentDist_ValuesOnGrid_InfHorz of that type (its Policy and its parameters), every entry, on a model with d and z, and on
% one where the number of z states differs by type (n_z a structure). (InfHorz counterpart of CoreFHorzPTypeTests parts 8 and 16.)

Names_i={'low','high'};
N_i=2;
n_d_withd=n_d; d_grid_withd=d_grid;
for cc=1:2
    FnsToEvaluate=struct(); vfoptions=struct();
    if cc==1
        mstr='with d and z';
        n_d=n_d_withd; d_grid=d_grid_withd;
        n_zs=[n_z,n_z]; n_zP=n_z; z_gridP=z_grid; pi_zP=pi_z;
        zg={z_grid,z_grid}; pz={pi_z,pi_z};
        ReturnFn=@(d,aprime,a,z,r,w_pt,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w_pt,sigma,eta,varphi);
        FnsToEvaluate.assets=@(d,aprime,a,z) a;
        FnsToEvaluate.earnings=@(d,aprime,a,z,w_pt) w_pt*d*z;
        FnsToEvaluate.savings=@(d,aprime,a,z,r) aprime-(1+r)*a;
    else
        mstr='n_z per type';
        n_d=0; d_grid=[];
        n_zs=[5,3];
        n_zP=struct(); z_gridP=struct(); pi_zP=struct(); zg=cell(N_i,1); pz=cell(N_i,1);
        for ii=1:N_i
            [zg{ii},pz{ii}]=discretizeAR1_FarmerToda(0,0.9,0.05,n_zs(ii)); zg{ii}=exp(zg{ii});
            n_zP.(Names_i{ii})=n_zs(ii); z_gridP.(Names_i{ii})=zg{ii}; pi_zP.(Names_i{ii})=pz{ii};
        end
        ReturnFn=@(aprime,a,z,r,w_pt,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w_pt,sigma);
        FnsToEvaluate.assets=@(aprime,a,z) a;
        FnsToEvaluate.earnings=@(aprime,a,z,w_pt) w_pt*z;
    end
    FnNames=fieldnames(FnsToEvaluate);
    [~,Policy_PT]=ValueFnIter_InfHorz_PType(n_d,n_a,n_zP,Names_i,d_grid,a_grid,z_gridP,pi_zP,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    VOG_PT=EvalFnOnAgentDist_ValuesOnGrid_InfHorz_PType(Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_zP,Names_i,d_grid,a_grid,z_gridP,struct());
    for ii=1:N_i
        nm=Names_i{ii};
        Params_ii=Params; Params_ii.w_pt=Params.w_pt(ii);
        VOG_ii=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy_PT.(nm),FnsToEvaluate,Params_ii,[],n_d,n_a,n_zs(ii),d_grid,a_grid,zg{ii},struct());
        for ff=1:numel(FnNames)
            x=gather(VOG_PT.(FnNames{ff}).(nm)); y=gather(VOG_ii.(FnNames{ff}));
            if ~isequal(size(x),size(y))
                fprintf('ValuesOnGrid PType %s, %s (type %s) size [%s] differs from single-type [%s], this should be zero: 1 \n',mstr,FnNames{ff},nm,num2str(size(x)),num2str(size(y)))
            else
                fprintf('ValuesOnGrid PType %s, %s (type %s) vs single-type, this should be zero: %.3e \n',mstr,FnNames{ff},nm,max(abs(x(:)-y(:))))
            end
        end
    end
end
output=struct();
end
