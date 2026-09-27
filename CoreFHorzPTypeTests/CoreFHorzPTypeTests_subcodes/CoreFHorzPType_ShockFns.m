function output=CoreFHorzPType_ShockFns(n_a,n_z,n_e,N_j,a_grid,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames)
% Test: shocks given by vfoptions.ExogShockFn and vfoptions.EiidShockFn under PType, two ways:
%   (1) one ExogShockFn/EiidShockFn for all types, whose parameters depend on ptype
%       (rho_z_pt, sigma_e_pt as length-N_i vectors)
%   (2) a struct of functions keyed by Names_i, each reading its own (type-free) parameters
% Both must give the same per-type V, Policy, StationaryDist and AllStats (every statistic) as solo
% solves given the grids explicitly (the grids being what the same functions return), and the
% grouped AllStats Mean must be the ptypeweights-weighted sum of the solo means.
% When ExogShockFn/EiidShockFn is used the raw z_grid, pi_z, e_grid and pi_e are ignored, so none
% are given here.

Names_i={'low','high'};
N_i=length(Names_i);

n_d=0;
d_grid=[];

Params.sigma_z=0.03;
Params.rho_z_pt=[0.9; 0.5];
Params.sigma_e_pt=[0.1; 0.2];
Params.rho_z_low=Params.rho_z_pt(1);
Params.rho_z_high=Params.rho_z_pt(2);
Params.sigma_e_low=Params.sigma_e_pt(1);
Params.sigma_e_high=Params.sigma_e_pt(2);

ReturnFn=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate.assets=@(aprime,a,z,e) a;
FnsToEvaluate.income=@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;
FnNames=fieldnames(FnsToEvaluate);

jequaloneDist=zeros(n_a,n_z,n_e,'gpuArray');
jequaloneDist(1,ceil(n_z/2),ceil(n_e/2))=1;

ptw=Params.(PTypeDistParamNames{1});

%% Solo solves, grids given explicitly
V_solo=struct(); Policy_solo=struct(); Dist_solo=struct(); AllStats_solo=struct();
for ii=1:N_i
    nm=Names_i{ii};
    [z_grid_ii,pi_z_ii]=CoreFHorzPType_ExogShockFn(Params.rho_z_pt(ii),Params.sigma_z);
    [e_grid_ii,pi_e_ii]=CoreFHorzPType_EiidShockFn(Params.sigma_e_pt(ii));
    vfopt_solo=struct();
    vfopt_solo.n_e=n_e;
    vfopt_solo.e_grid=e_grid_ii;
    vfopt_solo.pi_e=pi_e_ii;
    simopt_solo=vfopt_solo;
    [V_solo.(nm),Policy_solo.(nm)]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_ii,pi_z_ii,ReturnFn,Params,DiscountFactorParamNames,[],vfopt_solo);
    Dist_solo.(nm)=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_solo.(nm),n_d,n_a,n_z,N_j,pi_z_ii,Params,simopt_solo);
    AllStats_solo.(nm)=EvalFnOnAgentDist_AllStats_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_ii,simopt_solo);
end

%% The two PType cases
casenames={'one fn, per-type params','struct of fns'};
pairs=cell(0,3);
for cc=1:2
    vfopt_PT=struct();
    vfopt_PT.n_e=n_e;
    if cc==1
        vfopt_PT.ExogShockFn=@(rho_z_pt,sigma_z) CoreFHorzPType_ExogShockFn(rho_z_pt,sigma_z);
        vfopt_PT.EiidShockFn=@(sigma_e_pt) CoreFHorzPType_EiidShockFn(sigma_e_pt);
    elseif cc==2
        vfopt_PT.ExogShockFn.low =@(rho_z_low,sigma_z) CoreFHorzPType_ExogShockFn(rho_z_low,sigma_z);
        vfopt_PT.ExogShockFn.high=@(rho_z_high,sigma_z) CoreFHorzPType_ExogShockFn(rho_z_high,sigma_z);
        vfopt_PT.EiidShockFn.low =@(sigma_e_low) CoreFHorzPType_EiidShockFn(sigma_e_low);
        vfopt_PT.EiidShockFn.high=@(sigma_e_high) CoreFHorzPType_EiidShockFn(sigma_e_high);
    end
    simopt_PT=vfopt_PT;

    [V_PT,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,[],[],ReturnFn,Params,DiscountFactorParamNames,vfopt_PT);
    StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,Names_i,[],Params,simopt_PT);
    AllStats_PT=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,[],simopt_PT);

    for ii=1:N_i
        nm=Names_i{ii};
        fprintf('Shock fns (%s), V    (type %s), this should be zero: %.3e \n',casenames{cc},nm,max(abs(V_PT.(nm)(:)-V_solo.(nm)(:))))
        fprintf('Shock fns (%s), Pol  (type %s), this should be zero: %.3e \n',casenames{cc},nm,max(abs(Policy_PT.(nm)(:)-Policy_solo.(nm)(:))))
        fprintf('Shock fns (%s), Dist (type %s), this should be zero: %.3e \n',casenames{cc},nm,max(abs(StationaryDist_PT.(nm)(:)-Dist_solo.(nm)(:))))
        for ff=1:length(FnNames)
            pairs(end+1,:)={sprintf('Shock fns (%s), AllStats %s (type %s)',casenames{cc},FnNames{ff},nm), AllStats_PT.(FnNames{ff}).(nm), AllStats_solo.(nm).(FnNames{ff})}; %#ok<AGROW>
        end
    end
    for ff=1:length(FnNames)
        fn=FnNames{ff};
        agg_mean=ptw(1)*AllStats_solo.(Names_i{1}).(fn).Mean+ptw(2)*AllStats_solo.(Names_i{2}).(fn).Mean;
        fprintf('Shock fns (%s), AllStats %s grouped Mean, this should be zero: %.3e \n',casenames{cc},fn,abs(AllStats_PT.(fn).Mean-agg_mean))
    end
end

%% The two types really do differ
fprintf('Shock fns, V of the two types differ, this should NOT be zero: %.3e \n',max(abs(V_solo.(Names_i{1})(:)-V_solo.(Names_i{2})(:))))

% Compare every statistic of each (label, Test, Ref) pair. Ref supplies the list of stats; per-ptype
% substructures are skipped (they get pairs of their own). NaN is compared as a pattern, the count of
% entries that are NaN in one but not the other, because a NaN printed as the value of a check is
% not read by CoreSummary at all: the check would silently vanish from the count.
for pp=1:size(pairs,1)
    Test=pairs{pp,2};
    Ref=pairs{pp,3};
    statnames=fieldnames(Ref);
    for ss=1:length(statnames)
        sn=statnames{ss};
        if ischar(Ref.(sn)) || (isstruct(Ref.(sn)) && ~strcmp(sn,'MoreInequality'))
            continue
        end
        if strcmp(sn,'MoreInequality')
            subnames=fieldnames(Ref.MoreInequality);
            cmpnames=cell(length(subnames),1); xs=cmpnames; ys=cmpnames;
            for s2=1:length(subnames)
                cmpnames{s2}=['MoreInequality.',subnames{s2}];
                ys{s2}=Ref.MoreInequality.(subnames{s2});
                if isfield(Test,'MoreInequality') && isfield(Test.MoreInequality,subnames{s2})
                    xs{s2}=Test.MoreInequality.(subnames{s2});
                else
                    xs{s2}=[];
                end
            end
        else
            cmpnames={sn};
            ys={Ref.(sn)};
            if isfield(Test,sn)
                xs={Test.(sn)};
            else
                xs={[]};
            end
        end
        for s2=1:length(cmpnames)
            x=xs{s2}; y=ys{s2};
            if numel(x)~=numel(y)
                fprintf('%s, %s missing or wrongly sized, this should be zero: 1 \n',pairs{pp,1},cmpnames{s2})
            else
                x=gather(x(:)); y=gather(y(:));
                dd=abs(x-y); dd=dd(~isnan(dd));
                fprintf('%s, %s, this should be zero: %.3e, NaN-pattern mismatches %i \n',pairs{pp,1},cmpnames{s2},max([dd;0]),sum(isnan(x)~=isnan(y)))
            end
        end
    end
end

output=struct();

end
