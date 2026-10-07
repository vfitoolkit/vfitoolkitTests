function output=CoreInfHorzPType_SimPanel(n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT)
% Test: SimPanelValues_InfHorz_PType (new), three types, the middle one of zero weight, starting from each type's stationary
% distribution (InitialDist a structure with a field per type). The panel is random, so only:
%   exact: SimPanelValues.(fn) is simperiods-by-numbersims; no entry is NaN; the number of simulations of each type is
%     floor(ptweights*numbersims) (the remainder going to the first types, as in SimPanelValues_FHorz_Case1_PType), so the
%     zero-weight type has none (a FnsToEvaluate 'typeid' tells each simulation's type)
%   loose: each type's mean over the panel vs the PType AllStats of that type, and the panel mean vs the grouped AllStats
%   exact again with InitialDist one distribution for every type (a point mass), the panel of each type starting there
% (InfHorz counterpart of the SimPanel checks of CoreFHorzPTypeTests parts 8, 9 and 11.)

Names_i={'low','ghost','high'};
N_i=3;
ptw=[0.5;0;0.5];
Params.(PTypeDistParamNames{1})=ptw;
Params.w_pt=[1;1.2;1.4];
Params.typeid_pt=[1;2;3];
n_d=0; d_grid=[];
ReturnFn=@(aprime,a,z,r,w_pt,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w_pt,sigma);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w_pt) w_pt*z;
FnsToEvaluate.typeid=@(aprime,a,z,typeid_pt) typeid_pt;
[~,Policy_PT]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,struct());
Dist_PT=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,Names_i,pi_z,Params,simoptionsPT);
AllStats_PT=EvalFnOnAgentDist_AllStats_InfHorz_PType(Dist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,struct());

simoptions=struct(); simoptions.numbersims=10^4+1; simoptions.simperiods=100;
InitialDist=struct();
for ii=1:N_i
    InitialDist.(Names_i{ii})=Dist_PT.(Names_i{ii});
end
for form=1:2
    if form==2
        InitialDist=zeros(n_a,n_z,'gpuArray'); InitialDist(1,ceil(n_z/2))=1; % one distribution for every type
        fstr='one InitialDist';
    else
        fstr='InitialDist per type';
    end
    SP=SimPanelValues_InfHorz_PType(InitialDist,PTypeDistParamNames,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions);
    nsims=floor(ptw*simoptions.numbersims); extra=simoptions.numbersims-sum(nsims); nsims(1:extra)=nsims(1:extra)+1;
    bad=0;
    for fn=fieldnames(FnsToEvaluate)'
        bad=bad+~isequal(size(SP.(fn{1})),[simoptions.simperiods,simoptions.numbersims])+sum(isnan(gather(SP.(fn{1})(:))));
    end
    fprintf('SimPanel PType %s, sizes wrong + NaN entries, this should be zero: %i \n',fstr,bad)
    ids=gather(SP.typeid(1,:));
    cnt=[sum(ids==1);sum(ids==2);sum(ids==3)];
    fprintf('SimPanel PType %s, simulations per type minus floor(ptweights*numbersims) (remainder to the first types), this should be zero: %i \n',fstr,sum(abs(cnt-nsims)))
    fprintf('SimPanel PType %s, a simulation whose type changes over the panel, this should be zero: %i \n',fstr,sum(any(gather(SP.typeid)~=ids,1)))
    if form==1
        for ii=[1,3]
            pm=mean(reshape(gather(SP.assets(:,ids==ii)),[],1));
            fprintf('SimPanel PType %s, assets mean of type %s over the panel vs AllStats of the type, this should be close to zero: %.3e \n',fstr,Names_i{ii},abs(pm-gather(AllStats_PT.assets.(Names_i{ii}).Mean)))
        end
        fprintf('SimPanel PType %s, assets mean over the panel vs the grouped AllStats, this should be close to zero: %.3e \n',fstr,abs(mean(gather(SP.assets(:)))-gather(AllStats_PT.assets.Mean)))
    else
        a1=gather(SP.assets(1,:));
        fprintf('SimPanel PType %s, period-1 assets not the InitialDist point, this should be zero: %i \n',fstr,sum(a1~=a_grid(1)))
    end
end
output=struct();
end
