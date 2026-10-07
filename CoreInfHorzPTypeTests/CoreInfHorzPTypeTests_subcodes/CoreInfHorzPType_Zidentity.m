function output=CoreInfHorzPType_Zidentity(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT)
% Test: every grouped statistic of AllStats and AggVars PType, via z-identity. A PType model with N_i types is the same model as one
% without PType whose exogenous state is (z, type), with the identity as the transition of the type: the type's wage w_pt is the
% value of the second shock. Grouping across ptypes is then pooling across the second shock, so the grouped PType stats must equal
% the plain single-type stats of the (z, type) model, every statistic; and each type's stats must equal the single-type stats
% restricted to that type (type<ii>: the second shock equals w_pt(ii)). That is the exact reference for the grouped stats beyond
% the Mean. (InfHorz counterpart of CoreFHorzPTypeTests parts 2 and 10.)
% The (z, type) model is given the PType Policy and StationaryDist (stacked, the distribution of type ii times its ptweight), so
% the comparison is of the stats commands alone. Separately: solving the (z, type) model directly gives the same V (and the same
% Policy), and its stationary distribution, iterated from that stacked one (simoptions.initialdist: with the identity on the type
% the distribution does not mix across types, so the type masses must be set by the starting point), stays there.

Names_i={'low','high'};
N_i=2;
ptw=Params.(PTypeDistParamNames{1});
vfoptions=struct();

% (A) PType, with d and z
ReturnFn_A=@(d,aprime,a,z,r,w_pt,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w_pt,sigma,eta,varphi);
FnsToEvaluate_A.assets=@(d,aprime,a,z) a;
FnsToEvaluate_A.earnings=@(d,aprime,a,z,w_pt) w_pt*d*z;
FnsToEvaluate_A.savings=@(d,aprime,a,z) aprime-a;
FnNames=fieldnames(FnsToEvaluate_A);
[V_A,Policy_A]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_A,Params,DiscountFactorParamNames,vfoptions);
Dist_A=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy_A,n_d,n_a,n_z,Names_i,pi_z,Params,simoptionsPT);
AggVars_A=EvalFnOnAgentDist_AggVars_InfHorz_PType(Dist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,struct());
AllStats_A=EvalFnOnAgentDist_AllStats_InfHorz_PType(Dist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,struct());

% (B) No PType: the exogenous state is (z, type), z varying first; the identity on the type
n_z_B=[n_z,N_i];
z_grid_B=[z_grid; Params.w_pt];
pi_z_B=kron(eye(N_i),gather(pi_z));
ReturnFn_B=@(d,aprime,a,z,wt,r,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,wt,sigma,eta,varphi);
FnsToEvaluate_B.assets=@(d,aprime,a,z,wt) a;
FnsToEvaluate_B.earnings=@(d,aprime,a,z,wt) wt*d*z;
FnsToEvaluate_B.savings=@(d,aprime,a,z,wt) aprime-a;
simoptions_B=struct();
simoptions_B.conditionalrestrictions.type1=@(d,aprime,a,z,wt) (wt<1.2); % w_pt is 1 and 1.4
simoptions_B.conditionalrestrictions.type2=@(d,aprime,a,z,wt) (wt>1.2);
% the PType Policy and StationaryDist, stacked over the type
nP=size(Policy_A.(Names_i{1}),1);
Policy_B=zeros(nP,n_a,n_z,N_i,'gpuArray');
Dist_B=zeros(n_a,n_z,N_i,'gpuArray');
for ii=1:N_i
    Policy_B(:,:,:,ii)=reshape(Policy_A.(Names_i{ii}),[nP,n_a,n_z]);
    Dist_B(:,:,ii)=ptw(ii)*reshape(Dist_A.(Names_i{ii}),[n_a,n_z]);
end
AggVars_B=EvalFnOnAgentDist_AggVars_InfHorz(Dist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,d_grid,a_grid,z_grid_B,struct());
AllStats_B=EvalFnOnAgentDist_AllStats_InfHorz(Dist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,d_grid,a_grid,z_grid_B,simoptions_B);

% the (z, type) model solved directly
[V_Bs,Policy_Bs]=ValueFnIter_InfHorz(n_d,n_a,n_z_B,d_grid,a_grid,z_grid_B,pi_z_B,ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptions);
simoptions_Bs=simoptionsPT; simoptions_Bs.initialdist=Dist_B;
Dist_Bs=StationaryDist_InfHorz(Policy_Bs,n_d,n_a,n_z_B,pi_z_B,simoptions_Bs,Params,[]);
dV=0; nPolicy=0;
for ii=1:N_i
    dV=max(dV,max(abs(reshape(gather(V_Bs(:,:,ii)),[],1)-gather(V_A.(Names_i{ii})(:)))));
    nPolicy=nPolicy+sum(reshape(gather(Policy_Bs(:,:,:,ii)),[],1)~=gather(Policy_A.(Names_i{ii})(:)));
end
fprintf('z-identity stats, V of the (z, type) model solved directly vs PType, this should be close to zero: %.3e \n',dV)
fprintf('z-identity stats, Policy of the (z, type) model solved directly vs PType, entries that differ (information only): %i \n',nPolicy)
fprintf('z-identity stats, StationaryDist of the (z, type) model iterated from the stacked PType one, this should be close to zero: %.3e \n',max(abs(gather(Dist_Bs(:))-gather(Dist_B(:)))))

% Grouped and per type, every leaf
for ff=1:numel(FnNames)
    fn=FnNames{ff};
    fprintf('z-identity stats, AggVars %s grouped vs the (z, type) model, this should be zero: %.3e \n',fn,abs(gather(AggVars_A.(fn).Mean)-gather(AggVars_B.(fn).Mean)))
end
for tt=0:N_i
    maxd=0; mis=0;
    for ff=1:numel(FnNames)
        fn=FnNames{ff};
        if tt==0
            [p0,v0]=CoreStructLeaves(AllStats_B.(fn)); [p1,v1]=CoreStructLeaves(AllStats_A.(fn)); ustr='grouped';
            keep=cellfun(@(p) ~any(strcmp(strtok(p,'.'),Names_i)),p1); p1=p1(keep); v1=v1(keep); % (drop the per-type substructures of the PType output)
        else
            [p0,v0]=CoreStructLeaves(AllStats_B.(['type',num2str(tt)]).(fn)); [p1,v1]=CoreStructLeaves(AllStats_A.(fn).(Names_i{tt})); ustr=['type ',Names_i{tt}];
        end
        for k=1:numel(p0)
            [present,idx]=ismember(p0{k},p1);
            if present && isequal(size(v0{k}),size(v1{idx}))
                x=v0{k}; y=v1{idx}; mis=mis+sum(isnan(x(:))~=isnan(y(:))); d=abs(x(:)-y(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
            else
                mis=mis+1;
            end
        end
    end
    fprintf('z-identity stats, AllStats %s vs the (z, type) model, every leaf, this should be zero: %.3e, mismatches %i \n',ustr,maxd,mis)
end
for ii=1:N_i
    fprintf('z-identity stats, the type%i restriction RestrictedSampleMass vs ptweight, this should be zero: %.3e \n',ii,abs(gather(AllStats_B.(['type',num2str(ii)]).RestrictedSampleMass)-ptw(ii)))
end
output=struct();
end
