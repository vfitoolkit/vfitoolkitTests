function output=CoreFHorzPType_SimPanelIdentityz(n_a,N_j,a_grid,Params,DiscountFactorParamNames,PTypeDistParamNames,N_i)
% Test: SimPanelValues_FHorz_Case1_PType, exactly, with z.
% The with-z counterpart of CoreFHorzPType_SimPanelNoz. z has identity transitions, and each type
% starts from a point mass with type ii at z_grid(ii), so z never moves and every agent of a given
% type follows the same deterministic path. So the PType panel must be the per-type solo path
% repeated, with the types in order.
%
% jequaloneDist is input with ptype as a dimension, [n_a,n_z,N_i], and slice masses different from
% Params.(PTypeDistParamNames{1}). The slice masses override PTypeDistParamNames, so they set the
% number of columns of each type. The struct form of the same per-type distributions uses
% PTypeDistParamNames instead, and is checked alongside as the contrast.
%
% Types differ in the CRRA curvature (sigma_pt) and in their starting z. The FnsToEvaluate include
% the type's own sigma and z, so a column holding the wrong type is caught even where two paths coincide.

Params.sigma_pt=[2; 3];

n_d=0;
d_grid=[];
n_z=3; % not 2: with n_z=N_i a shared [n_z,n_z] pi_z is read as having ptype as its last dimension
z_grid=[0.8; 1; 1.2];
pi_z=eye(n_z); % identity: z never moves

vfoptions=struct();

ReturnFn_PT=@(aprime,a,z,r,w,kappa_j,sigma_pt,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma_pt,agej,Jr,pension);
ReturnFn_NoPT=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate_PT.assets=@(aprime,a,z) a;
FnsToEvaluate_PT.savings=@(aprime,a,z) aprime;
FnsToEvaluate_PT.shock=@(aprime,a,z) z;
FnsToEvaluate_PT.curvature=@(aprime,a,z,sigma_pt) sigma_pt;
FnsToEvaluate_NoPT.assets=@(aprime,a,z) a;
FnsToEvaluate_NoPT.savings=@(aprime,a,z) aprime;
FnsToEvaluate_NoPT.shock=@(aprime,a,z) z;
FnsToEvaluate_NoPT.curvature=@(aprime,a,z,sigma) sigma;
FnNames=fieldnames(FnsToEvaluate_PT);

% ptype as a dimension: type ii has no assets and starts at z_grid(ii), with mass slicemasses(ii)
slicemasses=[0.3; 0.7];
jequaloneDist_ptdim=zeros(n_a,n_z,N_i,'gpuArray');
for ii=1:N_i
    jequaloneDist_ptdim(1,ii,ii)=slicemasses(ii);
end

numbersims=1000;
simoptions=struct();
simoptions.numbersims=numbersims;

%% PType
[~,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions);
names_PT=fieldnames(Policy_PT);
SimPanel_ptdim=SimPanelValues_FHorz_Case1_PType(jequaloneDist_ptdim,PTypeDistParamNames,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,simoptions);
% the same per-type distributions as a struct (each of mass one), which uses PTypeDistParamNames
jequaloneDist_struct=struct();
for ii=1:N_i
    jequaloneDist_struct.(names_PT{ii})=jequaloneDist_ptdim(:,:,ii)/slicemasses(ii);
end
SimPanel_struct=SimPanelValues_FHorz_Case1_PType(jequaloneDist_struct,PTypeDistParamNames,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,simoptions);

%% Solo solves, a small panel each
simoptions_solo=struct();
simoptions_solo.numbersims=10;
SimPanel_solo=cell(N_i,1);
for ii=1:N_i
    Params_ii=Params;
    Params_ii.sigma=Params.sigma_pt(ii);
    jequaloneDist_ii=jequaloneDist_ptdim(:,:,ii)/slicemasses(ii);
    [~,Policy_ii]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_NoPT,Params_ii,DiscountFactorParamNames,[],vfoptions);
    SimPanel_solo{ii}=SimPanelValues_FHorz_Case1(jequaloneDist_ii,Policy_ii,FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_solo);
end

%% Expected columns of each type
% ptype-as-dim: from the slice masses
nsims_ptdim=floor(slicemasses*numbersims);
extrasims_ptdim=numbersims-sum(nsims_ptdim);
nsims_ptdim(1:extrasims_ptdim)=nsims_ptdim(1:extrasims_ptdim)+1;
lastcol_ptdim=cumsum(nsims_ptdim);
firstcol_ptdim=lastcol_ptdim-nsims_ptdim+1;
% struct: from PTypeDistParamNames
ptw=Params.(PTypeDistParamNames{1});
nsims_struct=floor(ptw*numbersims);
extrasims_struct=numbersims-sum(nsims_struct);
nsims_struct(1:extrasims_struct)=nsims_struct(1:extrasims_struct)+1;
lastcol_struct=cumsum(nsims_struct);
firstcol_struct=lastcol_struct-nsims_struct+1;

%% Compare
fprintf('SimPanel identity z, columns per type of ptype-as-dim differ from those of struct, this should NOT be zero: %i \n',max(abs(nsims_ptdim-nsims_struct)))
for ff=1:length(FnNames)
    fn=FnNames{ff};
    fprintf('SimPanel identity z (ptype-as-dim jequaloneDist), %s number of sims minus numbersims, this should be zero: %i \n',fn,size(SimPanel_ptdim.(fn),2)-numbersims)
    fprintf('SimPanel identity z (ptype-as-dim jequaloneDist), %s number of NaN entries, this should be zero: %i \n',fn,sum(isnan(SimPanel_ptdim.(fn)(:))))
    fprintf('SimPanel identity z (struct jequaloneDist), %s number of sims minus numbersims, this should be zero: %i \n',fn,size(SimPanel_struct.(fn),2)-numbersims)
    fprintf('SimPanel identity z (struct jequaloneDist), %s number of NaN entries, this should be zero: %i \n',fn,sum(isnan(SimPanel_struct.(fn)(:))))
    for ii=1:N_i
        solopath=gather(SimPanel_solo{ii}.(fn)(:,1));
        % The solo panel is itself deterministic (every column the same path)
        fprintf('SimPanel identity z, %s (noPType, type %i) columns all the same path, this should be zero: %.3e \n',fn,ii,max(max(abs(gather(SimPanel_solo{ii}.(fn))-solopath))))
        PTcols=gather(SimPanel_ptdim.(fn)(:,firstcol_ptdim(ii):lastcol_ptdim(ii)));
        fprintf('SimPanel identity z (ptype-as-dim jequaloneDist), %s (PType, %s) every column equals the solo path, this should be zero: %.3e \n',fn,names_PT{ii},max(max(abs(PTcols-solopath))))
        PTcols=gather(SimPanel_struct.(fn)(:,firstcol_struct(ii):lastcol_struct(ii)));
        fprintf('SimPanel identity z (struct jequaloneDist), %s (PType, %s) every column equals the solo path, this should be zero: %.3e \n',fn,names_PT{ii},max(max(abs(PTcols-solopath))))
    end
end

output=struct();

end
