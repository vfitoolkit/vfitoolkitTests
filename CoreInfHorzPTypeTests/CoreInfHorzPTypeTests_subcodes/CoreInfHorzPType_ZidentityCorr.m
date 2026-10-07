function output=CoreInfHorzPType_ZidentityCorr(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT)
% Test: the grouped outputs of the two new PType correlation commands
%   EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType
%   EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType
% via z-identity, as part 2: the PType model is the single-type model whose exogenous state is (z, type), the identity being the
% transition of the type, given the PType Policy and the stacked PType StationaryDist. Grouping across ptypes is pooling across
% the type, so (exactly, every field):
%   CrossSection grouped = single-type on the (z, type) model; per type = single-type restricted to the type (type<ii>); and the
%     same for the restriction rich (a>0.5): grouped = rich, per type = richtype<ii> (rich and of type ii)
%   AutoCorr grouped Mean/StdDeviation = single-type; restricted grouped (rich) = single-type rich, every field (pair fields
%     included: the grouped PairMass is a population share, as rich's); per type restricted = richtype<ii>, the PairMass rescaled
%     by the ptweight (the per-type PairMass is a within-type share)
% Up to the convergence error of the stationary distribution (solved to 1e-12 here, so still identities at the bar): the
% unrestricted grouped AutoCovariance/AutoCorrelation, as the PType grouping takes each type to be stationary while the
% single-type command pushes the stacked distribution forward.
% (InfHorz counterpart of CoreFHorzPTypeTests part 20.)

Names_i={'low','high'};
N_i=2;
ptw=Params.(PTypeDistParamNames{1});
timehorizons=[2,3]; horizons=[1,2,3];

% (A) PType, with d and z
ReturnFn_A=@(d,aprime,a,z,r,w_pt,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w_pt,sigma,eta,varphi);
FnsToEvaluate_A.assets=@(d,aprime,a,z) a;
FnsToEvaluate_A.earnings=@(d,aprime,a,z,w_pt) w_pt*d*z;
FnsToEvaluate_A.hours=@(d,aprime,a,z) d;
FnNames=fieldnames(FnsToEvaluate_A); nF=numel(FnNames);
[~,Policy_A]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_A,Params,DiscountFactorParamNames,struct());
Dist_A=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy_A,n_d,n_a,n_z,Names_i,pi_z,Params,simoptionsPT);
simoptions_A=struct(); simoptions_A.conditionalrestrictions.rich=@(d,aprime,a,z) (a>0.5); simoptions_A.timehorizons=timehorizons;
CS_A=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(Dist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_A);
ACP_A=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType(Dist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_A);

% (B) the (z, type) model, with the PType Policy and the stacked StationaryDist
n_z_B=[n_z,N_i];
z_grid_B=[z_grid; Params.w_pt];
pi_z_B=kron(eye(N_i),gather(pi_z));
FnsToEvaluate_B.assets=@(d,aprime,a,z,wt) a;
FnsToEvaluate_B.earnings=@(d,aprime,a,z,wt) wt*d*z;
FnsToEvaluate_B.hours=@(d,aprime,a,z,wt) d;
simoptions_B=struct(); simoptions_B.timehorizons=timehorizons;
simoptions_B.conditionalrestrictions.rich=@(d,aprime,a,z,wt) (a>0.5);
simoptions_B.conditionalrestrictions.type1=@(d,aprime,a,z,wt) (wt<1.2); % w_pt is 1 and 1.4
simoptions_B.conditionalrestrictions.type2=@(d,aprime,a,z,wt) (wt>1.2);
simoptions_B.conditionalrestrictions.richtype1=@(d,aprime,a,z,wt) (a>0.5)*(wt<1.2);
simoptions_B.conditionalrestrictions.richtype2=@(d,aprime,a,z,wt) (a>0.5)*(wt>1.2);
nP=size(Policy_A.(Names_i{1}),1);
Policy_B=zeros(nP,n_a,n_z,N_i,'gpuArray');
Dist_B=zeros(n_a,n_z,N_i,'gpuArray');
for ii=1:N_i
    Policy_B(:,:,:,ii)=reshape(Policy_A.(Names_i{ii}),[nP,n_a,n_z]);
    Dist_B(:,:,ii)=ptw(ii)*reshape(Dist_A.(Names_i{ii}),[n_a,n_z]);
end
CS_B=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(Dist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,d_grid,a_grid,z_grid_B,simoptions_B);
ACP_B=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(Dist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,d_grid,a_grid,z_grid_B,pi_z_B,simoptions_B);

%% CrossSection: grouped, per type, rich grouped, rich per type; every field of each (label, Test, Ref) pair
pairs=cell(0,3);
for ff=1:nF
    fn=FnNames{ff};
    pairs(end+1,:)={sprintf('z-identity corr, CrossSection %s grouped',fn), CS_A.(fn), CS_B.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity corr, CrossSection rich %s grouped',fn), CS_A.rich.(fn), CS_B.rich.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity corr, AutoCorr rich %s grouped',fn), ACP_A.rich.(fn), ACP_B.rich.(fn)}; %#ok<AGROW>
    for ii=1:N_i
        nm=Names_i{ii};
        pairs(end+1,:)={sprintf('z-identity corr, CrossSection %s (type %s)',fn,nm), CS_A.(fn).(nm), CS_B.(['type',num2str(ii)]).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity corr, CrossSection rich %s (type %s)',fn,nm), CS_A.rich.(fn).(nm), CS_B.(['richtype',num2str(ii)]).(fn)}; %#ok<AGROW>
        % the per-type PairMass is a within-type share, that of richtype<ii> a population share
        Ref=ACP_B.(['richtype',num2str(ii)]).(fn);
        Ref.PairMass=Ref.PairMass/ptw(ii);
        for hh=2:numel(horizons)
            Ref.(['tperiods',num2str(horizons(hh))]).PairMass=Ref.(['tperiods',num2str(horizons(hh))]).PairMass/ptw(ii);
        end
        pairs(end+1,:)={sprintf('z-identity corr, AutoCorr rich %s (type %s)',fn,nm), ACP_A.rich.(fn).(nm), Ref}; %#ok<AGROW>
    end
end
for pp=1:size(pairs,1)
    [p0,v0]=CoreStructLeaves(pairs{pp,3}); [p1,v1]=CoreStructLeaves(pairs{pp,2});
    maxd=0; mis=0;
    for k=1:numel(p0)
        [present,idx]=ismember(p0{k},p1);
        if present && isequal(size(v0{k}),size(v1{idx}))
            x=v0{k}; y=v1{idx}; mis=mis+sum(isnan(x(:))~=isnan(y(:))); d=abs(x(:)-y(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
        else
            mis=mis+1;
        end
    end
    fprintf('%s, every field of the reference, this should be zero: %.3e, mismatches %i \n',pairs{pp,1},maxd,mis)
end
x=[gather(CS_A.CovarianceMatrix(:));gather(CS_A.CorrelationMatrix(:));gather(CS_A.rich.CovarianceMatrix(:));gather(CS_A.rich.CorrelationMatrix(:))];
y=[gather(CS_B.CovarianceMatrix(:));gather(CS_B.CorrelationMatrix(:));gather(CS_B.rich.CovarianceMatrix(:));gather(CS_B.rich.CorrelationMatrix(:))];
d=abs(x-y); d=d(~isnan(d));
fprintf('z-identity corr, CrossSection grouped and rich grouped matrices, this should be zero: %.3e, NaN-pattern mismatches %i \n',max([d;0]),sum(isnan(x)~=isnan(y)))
fprintf('z-identity corr, CrossSection rich RestrictedSampleMass.Total vs the (z, type) model, this should be zero: %.3e \n',abs(gather(CS_A.rich.RestrictedSampleMass.Total)-gather(CS_B.rich.RestrictedSampleMass)))
fprintf('z-identity corr, AutoCorr rich RestrictedSampleMass.Total vs the (z, type) model, this should be zero: %.3e \n',abs(gather(ACP_A.rich.RestrictedSampleMass.Total)-gather(ACP_B.rich.RestrictedSampleMass)))

%% AutoCorr unrestricted grouped: Mean/StdDeviation exactly; AutoCovariance/AutoCorrelation to the convergence of the stationary dist
dm=0; da=0;
for ff=1:nF
    fn=FnNames{ff};
    dm=max([dm,abs(gather(ACP_A.(fn).Mean)-gather(ACP_B.(fn).Mean)),abs(gather(ACP_A.(fn).StdDeviation)-gather(ACP_B.(fn).StdDeviation))]);
    for hh=1:numel(horizons)
        if horizons(hh)==1
            GA=ACP_A.(fn); GB=ACP_B.(fn);
        else
            GA=ACP_A.(fn).(['tperiods',num2str(horizons(hh))]); GB=ACP_B.(fn).(['tperiods',num2str(horizons(hh))]);
        end
        da=max([da,abs(gather(GA.AutoCovariance)-gather(GB.AutoCovariance)),abs(gather(GA.AutoCorrelation)-gather(GB.AutoCorrelation))]);
    end
end
fprintf('z-identity corr, AutoCorr grouped Mean/StdDeviation vs the (z, type) model, this should be zero: %.3e \n',dm)
fprintf('z-identity corr, AutoCorr grouped AutoCovariance/AutoCorrelation (horizons 1,2,3) vs the (z, type) model (the stationary dist solved to 1e-12), this should be zero: %.3e \n',da)
output=struct();
end
