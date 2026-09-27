function output=CoreFHorzPType_ShockGridForms(n_a,n_z,n_e,N_j,a_grid,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames)
% Test: the two ways of making z_grid, pi_z, e_grid and pi_e depend on ptype.
%   struct          z_grid.(name), pi_z.(name), vfoptions.e_grid.(name), vfoptions.pi_e.(name)
%   trailing dim    the same arrays with ptype as an extra LAST dimension of length N_i
% each for age-independent and for age-dependent shocks, so four cases:
%   age-independent  z_grid [n_z,N_i]      pi_z [n_z,n_z,N_i]      e_grid [n_e,N_i]      pi_e [n_e,N_i]
%   age-dependent    z_grid [n_z,N_j,N_i]  pi_z [n_z,n_z,N_j,N_i]  e_grid [n_e,N_j,N_i]  pi_e [n_e,N_j,N_i]
% (struct forms: the same without the last dimension). The ShockTests (part 7) only ever use the
% struct form, and never age-dependent shocks.
% Each case must give the same per-type V, Policy, StationaryDist and AllStats (every statistic)
% as solo solves with that type's grids, and the grouped AllStats Mean must be the ptypeweights-
% weighted sum of the solo means.
%
% The two types have different z persistence and different e dispersion. In the age-dependent
% cases the grids also drift with age (differently for z and e), so reading a grid at the wrong
% age, or the wrong type's grid, changes the answer.

Names_i={'low','high'};
N_i=length(Names_i);

n_d=0;
d_grid=[];

rho_z_pt=[0.9; 0.5];
sigma_e_pt=[0.1; 0.2];

ReturnFn=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate.assets=@(aprime,a,z,e) a;
FnsToEvaluate.income=@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;
FnNames=fieldnames(FnsToEvaluate);

jequaloneDist=zeros(n_a,n_z,n_e,'gpuArray');
jequaloneDist(1,ceil(n_z/2),ceil(n_e/2))=1;

ptw=Params.(PTypeDistParamNames{1});

%% Per-type grids, age-independent and age-dependent
z_grid_I=struct(); pi_z_I=struct(); e_grid_I=struct(); pi_e_I=struct();
z_grid_J=struct(); pi_z_J=struct(); e_grid_J=struct(); pi_e_J=struct();
for ii=1:N_i
    nm=Names_i{ii};
    [zg,pz]=discretizeAR1_FarmerToda(0,rho_z_pt(ii),0.03,n_z);
    zg=exp(zg);
    [eg,pe]=discretizeAR1_FarmerToda(0,0,sigma_e_pt(ii),n_e);
    pe=pe(1,:)';
    eg=exp(eg);
    z_grid_I.(nm)=gpuArray(zg);
    pi_z_I.(nm)=gpuArray(pz);
    e_grid_I.(nm)=gpuArray(eg);
    pi_e_I.(nm)=gpuArray(pe);
    z_grid_J.(nm)=gpuArray(zg.*(1+0.02*(0:1:N_j-1)));
    pi_z_J.(nm)=gpuArray(repmat(pz,1,1,N_j));
    e_grid_J.(nm)=gpuArray(eg.*(1-0.01*(0:1:N_j-1)));
    pi_e_J.(nm)=gpuArray(repmat(pe,1,N_j));
end

%% Solo solves, age-independent (I) and age-dependent (J)
V_solo=struct(); Policy_solo=struct(); Dist_solo=struct(); AllStats_solo=struct();
agecases={'I','J'};
for cc=1:2
    ac=agecases{cc};
    for ii=1:N_i
        nm=Names_i{ii};
        if strcmp(ac,'I')
            z_grid_ii=z_grid_I.(nm); pi_z_ii=pi_z_I.(nm); e_grid_ii=e_grid_I.(nm); pi_e_ii=pi_e_I.(nm);
        else
            z_grid_ii=z_grid_J.(nm); pi_z_ii=pi_z_J.(nm); e_grid_ii=e_grid_J.(nm); pi_e_ii=pi_e_J.(nm);
        end
        vfopt_solo=struct();
        vfopt_solo.n_e=n_e;
        vfopt_solo.e_grid=e_grid_ii;
        vfopt_solo.pi_e=pi_e_ii;
        simopt_solo=vfopt_solo;
        [V_solo.(ac).(nm),Policy_solo.(ac).(nm)]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_ii,pi_z_ii,ReturnFn,Params,DiscountFactorParamNames,[],vfopt_solo);
        Dist_solo.(ac).(nm)=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_solo.(ac).(nm),n_d,n_a,n_z,N_j,pi_z_ii,Params,simopt_solo);
        AllStats_solo.(ac).(nm)=EvalFnOnAgentDist_AllStats_FHorz_Case1(Dist_solo.(ac).(nm),Policy_solo.(ac).(nm),FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_ii,simopt_solo);
    end
end

%% The four PType cases
casenames={'struct, age-independent','trailing dim, age-independent','struct, age-dependent','trailing dim, age-dependent'};
caseage={'I','I','J','J'};
pairs=cell(0,3);
for cc=1:4
    ac=caseage{cc};
    vfopt_PT=struct();
    vfopt_PT.n_e=n_e;
    if cc==1
        z_grid_PT=z_grid_I; pi_z_PT=pi_z_I;
        vfopt_PT.e_grid=e_grid_I; vfopt_PT.pi_e=pi_e_I;
    elseif cc==2
        z_grid_PT=cat(2,z_grid_I.(Names_i{1}),z_grid_I.(Names_i{2}));   % [n_z,N_i]
        pi_z_PT=cat(3,pi_z_I.(Names_i{1}),pi_z_I.(Names_i{2}));         % [n_z,n_z,N_i]
        vfopt_PT.e_grid=cat(2,e_grid_I.(Names_i{1}),e_grid_I.(Names_i{2})); % [n_e,N_i]
        vfopt_PT.pi_e=cat(2,pi_e_I.(Names_i{1}),pi_e_I.(Names_i{2}));       % [n_e,N_i]
    elseif cc==3
        z_grid_PT=z_grid_J; pi_z_PT=pi_z_J;
        vfopt_PT.e_grid=e_grid_J; vfopt_PT.pi_e=pi_e_J;
    elseif cc==4
        z_grid_PT=cat(3,z_grid_J.(Names_i{1}),z_grid_J.(Names_i{2}));   % [n_z,N_j,N_i]
        pi_z_PT=cat(4,pi_z_J.(Names_i{1}),pi_z_J.(Names_i{2}));         % [n_z,n_z,N_j,N_i]
        vfopt_PT.e_grid=cat(3,e_grid_J.(Names_i{1}),e_grid_J.(Names_i{2})); % [n_e,N_j,N_i]
        vfopt_PT.pi_e=cat(3,pi_e_J.(Names_i{1}),pi_e_J.(Names_i{2}));       % [n_e,N_j,N_i]
    end
    simopt_PT=vfopt_PT;

    [V_PT,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid_PT,pi_z_PT,ReturnFn,Params,DiscountFactorParamNames,vfopt_PT);
    StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,Names_i,pi_z_PT,Params,simopt_PT);
    AllStats_PT=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid_PT,simopt_PT);

    for ii=1:N_i
        nm=Names_i{ii};
        fprintf('Shock grid forms (%s), V    (type %s), this should be zero: %.3e \n',casenames{cc},nm,max(abs(V_PT.(nm)(:)-V_solo.(ac).(nm)(:))))
        fprintf('Shock grid forms (%s), Pol  (type %s), this should be zero: %.3e \n',casenames{cc},nm,max(abs(Policy_PT.(nm)(:)-Policy_solo.(ac).(nm)(:))))
        fprintf('Shock grid forms (%s), Dist (type %s), this should be zero: %.3e \n',casenames{cc},nm,max(abs(StationaryDist_PT.(nm)(:)-Dist_solo.(ac).(nm)(:))))
        for ff=1:length(FnNames)
            pairs(end+1,:)={sprintf('Shock grid forms (%s), AllStats %s (type %s)',casenames{cc},FnNames{ff},nm), AllStats_PT.(FnNames{ff}).(nm), AllStats_solo.(ac).(nm).(FnNames{ff})}; %#ok<AGROW>
        end
    end
    for ff=1:length(FnNames)
        fn=FnNames{ff};
        agg_mean=ptw(1)*AllStats_solo.(ac).(Names_i{1}).(fn).Mean+ptw(2)*AllStats_solo.(ac).(Names_i{2}).(fn).Mean;
        fprintf('Shock grid forms (%s), AllStats %s grouped Mean, this should be zero: %.3e \n',casenames{cc},fn,abs(AllStats_PT.(fn).Mean-agg_mean))
    end
end

%% The two types really do differ (else the per-type checks cannot tell the types' grids apart)
fprintf('Shock grid forms, V of the two types differ (age-independent), this should NOT be zero: %.3e \n',max(abs(V_solo.I.(Names_i{1})(:)-V_solo.I.(Names_i{2})(:))))
fprintf('Shock grid forms, V age-dependent differs from age-independent (type %s), this should NOT be zero: %.3e \n',Names_i{1},max(abs(V_solo.I.(Names_i{1})(:)-V_solo.J.(Names_i{1})(:))))

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
