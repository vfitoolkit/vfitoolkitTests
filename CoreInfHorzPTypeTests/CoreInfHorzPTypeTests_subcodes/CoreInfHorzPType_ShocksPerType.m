function output=CoreInfHorzPType_ShocksPerType(n_a,a_grid,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT)
% Test: the exogenous shock differs by type, four ways, each against separate single-type solves of each type (V, Policy,
% StationaryDist, AllStats every leaf, AggVars), and the grouped AllStats Mean/StdDeviation against the ptweights-pooling:
%   (a) z_grid and pi_z as structures with a field per type (the same n_z)
%   (b) n_z, z_grid and pi_z as structures with a field per type, a different number of states per type
%   (c) one vfoptions.ExogShockFn for all types, whose parameter (rho_z_pt) differs by type
%   (d) vfoptions.ExogShockFn as a structure with a function per type
% (InfHorz counterpart of CoreFHorzPTypeTests parts 13, 14 and 16.)

Names_i={'low','high'};
N_i=2;
ptw=Params.(PTypeDistParamNames{1});
n_d=0; d_grid=[];
ReturnFn=@(aprime,a,z,r,w_pt,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w_pt,sigma);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w_pt) w_pt*z;
FnNames=fieldnames(FnsToEvaluate);
Params.sigma_z=0.05;
Params.rho_z_pt=[0.8;0.95];
Params.rho_z_low=0.8; Params.rho_z_high=0.95;
n_zs=[5,5]; % the number of z states of each type (5 is what CoreInfHorzPType_ExogShockFn gives)

for form=1:4
    vfoptions=struct(); simoptions=simoptionsPT;
    if form==1
        fstr='(a) z_grid and pi_z per type';
        n_zs=[5,5];
        n_z=5; z_grid=struct(); pi_z=struct();
    elseif form==2
        fstr='(b) n_z, z_grid and pi_z per type, different sizes';
        n_zs=[5,3];
        n_z=struct(); z_grid=struct(); pi_z=struct();
    elseif form==3
        fstr='(c) one ExogShockFn, a per-type parameter';
        n_zs=[5,5];
        n_z=5; z_grid=[]; pi_z=[];
        vfoptions.ExogShockFn=@(rho_z_pt,sigma_z) CoreInfHorzPType_ExogShockFn(rho_z_pt,sigma_z);
        simoptions.ExogShockFn=vfoptions.ExogShockFn;
    else
        fstr='(d) ExogShockFn per type';
        n_zs=[5,5];
        n_z=5; z_grid=[]; pi_z=[];
        vfoptions.ExogShockFn.low=@(rho_z_low,sigma_z) CoreInfHorzPType_ExogShockFn(rho_z_low,sigma_z);
        vfoptions.ExogShockFn.high=@(rho_z_high,sigma_z) CoreInfHorzPType_ExogShockFn(rho_z_high,sigma_z);
        simoptions.ExogShockFn=vfoptions.ExogShockFn;
    end
    zg=cell(N_i,1); pz=cell(N_i,1);
    for ii=1:N_i
        [zg{ii},pz{ii}]=discretizeAR1_FarmerToda(0,Params.rho_z_pt(ii),Params.sigma_z,n_zs(ii));
        zg{ii}=exp(zg{ii});
        if form<=2
            z_grid.(Names_i{ii})=zg{ii}; pi_z.(Names_i{ii})=pz{ii};
            if form==2
                n_z.(Names_i{ii})=n_zs(ii);
            end
        end
    end
    [V_PT,Policy_PT]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    Dist_PT=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,Names_i,pi_z,Params,simoptions);
    AggVars_PT=EvalFnOnAgentDist_AggVars_InfHorz_PType(Dist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions);
    AllStats_PT=EvalFnOnAgentDist_AllStats_InfHorz_PType(Dist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions);
    m=zeros(numel(FnNames),N_i); s=zeros(numel(FnNames),N_i);
    for ii=1:N_i
        nm=Names_i{ii};
        Params_ii=Params; Params_ii.w_pt=Params.w_pt(ii);
        [V_ii,Policy_ii]=ValueFnIter_InfHorz(n_d,n_a,n_zs(ii),d_grid,a_grid,zg{ii},pz{ii},ReturnFn,Params_ii,DiscountFactorParamNames,[],struct());
        Dist_ii=StationaryDist_InfHorz(Policy_ii,n_d,n_a,n_zs(ii),pz{ii},simoptionsPT,Params_ii,[]);
        AggVars_ii=EvalFnOnAgentDist_AggVars_InfHorz(Dist_ii,Policy_ii,FnsToEvaluate,Params_ii,[],n_d,n_a,n_zs(ii),d_grid,a_grid,zg{ii},struct());
        AllStats_ii=EvalFnOnAgentDist_AllStats_InfHorz(Dist_ii,Policy_ii,FnsToEvaluate,Params_ii,[],n_d,n_a,n_zs(ii),d_grid,a_grid,zg{ii},struct());
        fprintf('Shocks per type %s, V (type %s) vs single-type, this should be zero: %.3e \n',fstr,nm,max(abs(gather(V_PT.(nm)(:))-gather(V_ii(:)))))
        fprintf('Shocks per type %s, Policy (type %s) vs single-type, this should be zero: %.3e \n',fstr,nm,max(abs(gather(Policy_PT.(nm)(:))-gather(Policy_ii(:)))))
        fprintf('Shocks per type %s, StationaryDist (type %s) vs single-type, this should be zero: %.3e \n',fstr,nm,max(abs(gather(Dist_PT.(nm)(:))-gather(Dist_ii(:)))))
        maxd=0; mis=0;
        for ff=1:numel(FnNames)
            fn=FnNames{ff};
            fprintf('Shocks per type %s, AggVars %s (type %s) vs single-type, this should be zero: %.3e \n',fstr,fn,nm,abs(gather(AggVars_PT.(fn).(nm).Mean)-gather(AggVars_ii.(fn).Mean)))
            [p0,v0]=CoreStructLeaves(AllStats_ii.(fn)); [p1,v1]=CoreStructLeaves(AllStats_PT.(fn).(nm));
            for k=1:numel(p0)
                [present,idx]=ismember(p0{k},p1);
                if present && isequal(size(v0{k}),size(v1{idx}))
                    x=v0{k}; y=v1{idx}; mis=mis+sum(isnan(x(:))~=isnan(y(:))); d=abs(x(:)-y(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                else
                    mis=mis+1;
                end
            end
            m(ff,ii)=gather(AllStats_ii.(fn).Mean); s(ff,ii)=gather(AllStats_ii.(fn).StdDeviation);
        end
        fprintf('Shocks per type %s, AllStats (type %s) vs single-type, every leaf, this should be zero: %.3e, mismatches %i \n',fstr,nm,maxd,mis)
    end
    for ff=1:numel(FnNames)
        mG=sum(ptw'.*m(ff,:)); sG=sqrt(sum(ptw'.*(s(ff,:).^2+(m(ff,:)-mG).^2)));
        fprintf('Shocks per type %s, AllStats %s grouped Mean/StdDeviation vs pooled single-type, this should be zero: %.3e \n',fstr,FnNames{ff},max(abs(gather(AllStats_PT.(FnNames{ff}).Mean)-mG),abs(gather(AllStats_PT.(FnNames{ff}).StdDeviation)-sG)))
    end
end
output=struct();
end
