function output=CoreFHorzPType_SimPanelNoz(n_a,N_j,a_grid,Params,DiscountFactorParamNames,PTypeDistParamNames,N_i)
% Test: SimPanelValues_FHorz_Case1_PType, exactly.
% A panel is random, so it cannot normally be compared exactly (and SimPanel is not even
% reproducible run to run: parfor workers ignore the client rng). But with no exogenous shocks
% and a jequaloneDist that is a point mass, every agent of a given type follows the same
% deterministic path. So the PType panel must be the per-type solo path repeated, with the
% types in order: the first PType_numbersims(1) columns are type 1, the next are type 2.
% The number of columns per type is floor(ptypeweight*numbersims), with the leftover sims added
% one each to the first types; that rule is replicated here.
%
% Types differ in the CRRA curvature (sigma_pt), as in part 2. The FnsToEvaluate include the
% type's own sigma, so a column holding the wrong type is caught even where two paths coincide.

Params.sigma_pt=[2; 3];

n_d=0;
d_grid=[];
n_z=0;
z_grid=[];
pi_z=[];

vfoptions=struct();

ReturnFn_PT=@(aprime,a,r,w,kappa_j,sigma_pt,agej,Jr,pension) ...
    ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma_pt,agej,Jr,pension);
ReturnFn_NoPT=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate_PT.assets=@(aprime,a) a;
FnsToEvaluate_PT.savings=@(aprime,a) aprime;
FnsToEvaluate_PT.curvature=@(aprime,a,sigma_pt) sigma_pt;
FnsToEvaluate_NoPT.assets=@(aprime,a) a;
FnsToEvaluate_NoPT.savings=@(aprime,a) aprime;
FnsToEvaluate_NoPT.curvature=@(aprime,a,sigma) sigma;
FnNames=fieldnames(FnsToEvaluate_PT);

jequaloneDist=zeros(n_a,1,'gpuArray');
jequaloneDist(1)=1; % no assets

numbersims=1000;
simoptions=struct();
simoptions.numbersims=numbersims;

%% PType
[~,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions);
SimPanel_PT=SimPanelValues_FHorz_Case1_PType(jequaloneDist,PTypeDistParamNames,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,simoptions);
names_PT=fieldnames(Policy_PT);

%% Solo solves, a small panel each
simoptions_solo=struct();
simoptions_solo.numbersims=10;
SimPanel_solo=cell(N_i,1);
for ii=1:N_i
    Params_ii=Params;
    Params_ii.sigma=Params.sigma_pt(ii);
    [~,Policy_ii]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_NoPT,Params_ii,DiscountFactorParamNames,[],vfoptions);
    SimPanel_solo{ii}=SimPanelValues_FHorz_Case1(jequaloneDist,Policy_ii,FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_solo);
end

%% Expected columns of each type
ptw=Params.(PTypeDistParamNames{1});
nsims_ii=floor(ptw*numbersims);
extrasims=numbersims-sum(nsims_ii);
nsims_ii(1:extrasims)=nsims_ii(1:extrasims)+1;
lastcol=cumsum(nsims_ii);
firstcol=lastcol-nsims_ii+1;

%% Compare
for ff=1:length(FnNames)
    fn=FnNames{ff};
    fprintf('SimPanel no z, %s number of sims minus numbersims, this should be zero: %i \n',fn,size(SimPanel_PT.(fn),2)-numbersims)
    fprintf('SimPanel no z, %s number of NaN entries, this should be zero: %i \n',fn,sum(isnan(SimPanel_PT.(fn)(:))))
    for ii=1:N_i
        solopath=gather(SimPanel_solo{ii}.(fn)(:,1));
        % The solo panel is itself deterministic (every column the same path)
        fprintf('SimPanel no z, %s (noPType, type %i) columns all the same path, this should be zero: %.3e \n',fn,ii,max(max(abs(gather(SimPanel_solo{ii}.(fn))-solopath))))
        PTcols=gather(SimPanel_PT.(fn)(:,firstcol(ii):lastcol(ii)));
        fprintf('SimPanel no z, %s (PType, %s) every column equals the solo path, this should be zero: %.3e \n',fn,names_PT{ii},max(max(abs(PTcols-solopath))))
    end
end

output=struct();

end
