function output=CoreFHorzPType_AutoCorrBruteForce(n_a,n_z,n_e,N_j,a_grid,z_grid,pi_z,e_grid,pi_e,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames)
% Test: EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType (added 2026-09-29) against a by-hand pooling of brute-force
% joint distributions. (The single-type command is tested against the same brute force in CoreFHorzTests part 12,
% TestFnsToEvaluate_AutoCorrBruteForce and TestFnsToEvaluate_AutoCorrBruteForce_semiz.)
%
% Brute force: for each type and age jj the one-period transition kernel T_jj (N_states x N_states, dense) is built
% by loops from the Policy indexes and pi_z/pi_e, the k-period kernel is the product T_j0*...*T_{j0+k-1}, and the
% (unnormalised) joint distribution of (state at j0, state at j0+k) is dist_j0(s)*Pk(s,s'), weighted by the type's
% ptweight and pooled over the types; every grouped moment is then a plain sum over that pooled joint distribution,
% and under a restriction over its pairs (the restriction satisfied at both ages).
%
% Model: no d, markov z, iid e (state (a,z,e), N_states=n_a*n_z*n_e); two types with different kappa_j and DIFFERENT
% age weights (so the pooling weights ptweights*(type's mass at j0) matter), ptweights 0.6/0.4. Checks, all exact:
%   per-type outputs (every field, unrestricted and restricted) equal the solo command, with groupptypesforstats=0 too;
%   grouped Mean/StdDeviation (unrestricted and restricted) and grouped AutoCovariance/AutoCorrelation and the pair
%   fields at horizons 1,2,3 equal the by-hand pooling; RestrictedSampleMass.ByAge/.Total by hand;
%   a type of zero weight leaves the grouped output equal to the other type alone;
%   groupptypesforstats=0 drops the grouped fields;
%   per-type TransitionProbs (assets, earnings; requested through simoptions.transprobs) equal the solo command, and
%   those of assets the binned dense kernel of that type (rows summing to one where the origin bin has mass), with
%   their TransitionValues_j/_jplus1 labels and TransitionMass_j;
%   grouped TransitionProbs (assets, earnings; unique values with labels and masses, and quintiles) equal a by-hand
%   pooling of the types' dense-kernel joints with weights ptweights*(type's mass at the age), over the union of the
%   types' value bins, the quintile bins being those of the pooled age-j and age-(j+1) value distributions.
% Zero checks are printed as %.3e; NaN is compared as a pattern (count of entries NaN in one but not the other).

fprintf('\n========== CoreFHorzPType_AutoCorrBruteForce ==========\n')
n_d=0;
d_grid=[];
N_s=n_a*n_z*n_e;

vfoptions=struct();
vfoptions.n_e=n_e; vfoptions.e_grid=e_grid; vfoptions.pi_e=pi_e;
simoptions=struct();
simoptions.n_e=n_e; simoptions.e_grid=e_grid; simoptions.pi_e=pi_e;

ReturnFn=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate.assets=@(aprime,a,z,e) a;
FnsToEvaluate.earnings=@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;
FnsToEvaluate.consumption=@(aprime,a,z,e,r,w,kappa_j) (1+r)*a+w*kappa_j*z*e-aprime;
FnsToEvaluate.evar=@(aprime,a,z,e) e;
FnNames=fieldnames(FnsToEvaluate);

Restrictions.rich=@(aprime,a,z,e) (a>0.5);
Restrictions.highz=@(aprime,a,z,e) (z>1);
RestNames=fieldnames(Restrictions);

jequaloneDist=zeros(n_a,n_z,n_e,'gpuArray');
jequaloneDist(1,ceil(n_z/2),ceil(n_e/2))=1;


pi_z_cpu=gather(pi_z);
pi_e_cpu=gather(pi_e(:));

%% B. PType: two types with different kappa_j and different age weights
Names_i={'low','high'};
N_i=2;
Params_PT=Params;
Params_PT.kappa_j_pt=[Params.kappa_j; 1.3*Params.kappa_j];
Params_PT.mewj=struct();
Params_PT.mewj.low=Params.mewj;
Params_PT.mewj.high=0.9.^(0:N_j-1); Params_PT.mewj.high=Params_PT.mewj.high/sum(Params_PT.mewj.high);
Params_PT.(PTypeDistParamNames{1})=[0.6;0.4];
ptw=[0.6;0.4];
ReturnFn_PT=@(aprime,a,z,e,r,w,kappa_j_pt,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j_pt,sigma,agej,Jr,pension);
FnsToEvaluate_PT.assets=@(aprime,a,z,e) a;
FnsToEvaluate_PT.earnings=@(aprime,a,z,e,w,kappa_j_pt) w*kappa_j_pt*z*e;
FnsToEvaluate_PT.consumption=@(aprime,a,z,e,r,w,kappa_j_pt) (1+r)*a+w*kappa_j_pt*z*e-aprime;
FnsToEvaluate_PT.evar=@(aprime,a,z,e) e;
timehorizons_PT=[2,3];
horizons_PT=[1,2,3];
hstr_PT={'','_k2','_k3'};

simoptions_PT=simoptions;
simoptions_PT.timehorizons=timehorizons_PT;
simoptions_PT.conditionalrestrictions=Restrictions;
simoptions_PT.transprobs={'assets','earnings'}; % earnings differs by type (kappa_j_pt), so the union of the types' value bins in the grouped TransitionProbs is non-trivial
[~,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params_PT,DiscountFactorParamNames,vfoptions);
StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_PT,simoptions);
ACP_PT=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_PT);
simoptions_PT_ng=simoptions_PT; simoptions_PT_ng.groupptypesforstats=0;
ACP_PT_ng=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_PT_ng);

% Solo solves, and their brute-force joint distributions
ACP_solo=struct(); Dist_solo=struct(); Vals_solo=struct(); Rvals_solo=struct(); T_solo=struct();
for ii=1:N_i
    nm=Names_i{ii};
    Params_ii=Params;
    Params_ii.kappa_j=Params_PT.kappa_j_pt(ii,:);
    Params_ii.mewj=Params_PT.mewj.(nm);
    [~,Policy_ii]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
    Dist_solo.(nm)=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_ii,n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions);
    ACP_solo.(nm)=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(Dist_solo.(nm),Policy_ii,FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_PT);
    fprintf('AutoCorr brute force B, Dist (type %s) PType vs solo, this should be zero: %.3e \n',nm,max(abs(gather(StationaryDist_PT.(nm)(:))-gather(Dist_solo.(nm)(:)))))
    ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy_ii,FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    Vals_solo.(nm)=cell(length(FnNames),1);
    for ff=1:length(FnNames)
        Vals_solo.(nm){ff}=gpuArray(reshape(gather(ValuesOnGrid.(FnNames{ff})),[N_s,N_j]));
    end
    RestOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy_ii,Restrictions,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    Rvals_solo.(nm)=cell(length(RestNames),1);
    for rr=1:length(RestNames)
        Rvals_solo.(nm){rr}=gpuArray(double(reshape(gather(RestOnGrid.(RestNames{rr})),[N_s,N_j])~=0));
    end
    Policy_cpu=gather(Policy_ii);
    T_solo.(nm)=cell(N_j-1,1);
    for jj=1:N_j-1
        Tjj=zeros(N_s,N_s);
        for zz=1:n_z
            for ee=1:n_e
                rows=(1:n_a)'+n_a*(zz-1)+n_a*n_z*(ee-1);
                aplow=reshape(Policy_cpu(1,:,zz,ee,jj),[n_a,1]);
                for zz2=1:n_z
                    for ee2=1:n_e
                        cols=aplow+n_a*(zz2-1)+n_a*n_z*(ee2-1);
                        Tjj(rows+(cols-1)*N_s)=Tjj(rows+(cols-1)*N_s)+pi_z_cpu(zz,zz2)*pi_e_cpu(ee2);
                    end
                end
            end
        end
        T_solo.(nm){jj}=gpuArray(Tjj);
    end
end

% Per type vs solo, every field (unrestricted and restricted)
for ii=1:N_i
    nm=Names_i{ii};
    for ff=1:length(FnNames)
        fn=FnNames{ff};
        statnames=fieldnames(ACP_solo.(nm).(fn));
        for ss=1:length(statnames)
            if iscell(ACP_solo.(nm).(fn).(statnames{ss})) % TransitionProbs and its labels/masses: cells over ages, compared entry by entry
                yc=ACP_solo.(nm).(fn).(statnames{ss});
                for tt=1:2
                    if tt==1
                        xc=ACP_PT.(fn).(nm).(statnames{ss}); lbl='PType vs solo';
                    else
                        xc=ACP_PT_ng.(fn).(nm).(statnames{ss}); lbl='groupptypesforstats=0 vs solo';
                    end
                    dcell=0; nbad=0;
                    for cc=1:numel(yc)
                        if numel(xc{cc})~=numel(yc{cc})
                            nbad=nbad+1;
                        elseif ~isempty(yc{cc})
                            dcell=max(dcell,max(abs(gather(xc{cc}(:))-gather(yc{cc}(:)))));
                        end
                    end
                    fprintf('AutoCorr brute force B, %s %s (type %s) %s, this should be zero: %.3e, wrongly sized cells %i \n',fn,statnames{ss},nm,lbl,dcell,nbad)
                end
            else
                x=gather(ACP_PT.(fn).(nm).(statnames{ss})(:)); y=gather(ACP_solo.(nm).(fn).(statnames{ss})(:)); dd=abs(x-y); dd=dd(~isnan(dd));
                fprintf('AutoCorr brute force B, %s %s (type %s) PType vs solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,statnames{ss},nm,max([dd;0]),sum(isnan(x)~=isnan(y)))
                x=gather(ACP_PT_ng.(fn).(nm).(statnames{ss})(:)); dd=abs(x-y); dd=dd(~isnan(dd));
                fprintf('AutoCorr brute force B, %s %s (type %s) groupptypesforstats=0 vs solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,statnames{ss},nm,max([dd;0]),sum(isnan(x)~=isnan(y)))
            end
        end
        for rr=1:length(RestNames)
            rn=RestNames{rr};
            statnames=fieldnames(ACP_solo.(nm).(rn).(fn));
            for ss=1:length(statnames)
                x=gather(ACP_PT.(rn).(fn).(nm).(statnames{ss})(:)); y=gather(ACP_solo.(nm).(rn).(fn).(statnames{ss})(:)); dd=abs(x-y); dd=dd(~isnan(dd));
                fprintf('AutoCorr brute force B, %s %s %s (type %s) PType vs solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',rn,fn,statnames{ss},nm,max([dd;0]),sum(isnan(x)~=isnan(y)))
            end
        end
    end
    for rr=1:length(RestNames)
        rn=RestNames{rr};
        fprintf('AutoCorr brute force B, %s RestrictedSampleMass (type %s) PType vs solo, this should be zero: %.3e \n',rn,nm,max(abs(gather(ACP_PT.(rn).RestrictedSampleMass.(nm)(:))-gather(ACP_solo.(nm).(rn).RestrictedSampleMass(:)))))
    end
end
for ff=1:length(FnNames)
    fprintf('AutoCorr brute force B, %s groupptypesforstats=0 has no grouped Mean, this should be zero: %i \n',FnNames{ff},isfield(ACP_PT_ng.(FnNames{ff}),'Mean'))
end

% Per-type TransitionProbs (assets) vs the binned dense kernel of that type: row b of TransitionProbs{jj} is the age-jj
% joint distribution dist_jj(s)*T_jj(s,s') summed over the states in origin bin b and in destination bin b', divided by
% the origin bin's mass (the bins are the unique asset values at each age, which TransitionValues_j/_jplus1 label, and
% TransitionMass_j is the origin bins' within-age mass). The command puts a zero row where an origin bin has no mass,
% and normalises the row of every bin with positive mass (however small; GI leaves bins with mass below eps).
ffa=find(strcmp(FnNames,'assets'));
for ii=1:N_i
    nm=Names_i{ii};
    Dist_u_ii=gpuArray(reshape(gather(Dist_solo.(nm)),[N_s,N_j]));
    dtp=0; drow=0; dlab=0; dmass=0; nanmm=0;
    for jj=1:N_j-1
        massj=gather(sum(Dist_u_ii(:,jj)));
        if massj==0
            continue
        end
        x=gather(Vals_solo.(nm){ffa}(:,jj)); y=gather(Vals_solo.(nm){ffa}(:,jj+1));
        [vx,~,ix]=unique(x); [vy,~,iy]=unique(y);
        Sx=gpuArray(full(sparse(1:N_s,ix,1,N_s,length(vx)))); Sy=gpuArray(full(sparse(1:N_s,iy,1,N_s,length(vy))));
        distj=Dist_u_ii(:,jj)/massj;
        massbin=Sx'*distj;
        massbin_safe=massbin; massbin_safe(massbin==0)=1; % the command puts a zero row where the origin bin has no mass, and normalises every bin with positive mass
        P_bf=gather((Sx'*(distj.*T_solo.(nm){jj})*Sy)./massbin_safe);
        P=gather(ACP_PT.assets.(nm).TransitionProbs{jj});
        if ~isequal(size(P),size(P_bf))
            dtp=max(dtp,1);
        else
            dd=abs(P-P_bf); dd=dd(~isnan(dd)); dtp=max(dtp,max([dd(:);0])); nanmm=nanmm+sum(isnan(P(:))~=isnan(P_bf(:)));
            rs=sum(P,2); drow=max(drow,max(abs(rs(gather(massbin)>0)-1)));
        end
        lx=gather(ACP_PT.assets.(nm).TransitionValues_j{jj}(:)); ly=gather(ACP_PT.assets.(nm).TransitionValues_jplus1{jj}(:)); lm=gather(ACP_PT.assets.(nm).TransitionMass_j{jj}(:));
        if numel(lx)~=numel(vx) || numel(ly)~=numel(vy) || numel(lm)~=numel(vx)
            dlab=max(dlab,1); dmass=max(dmass,1);
        else
            dlab=max([dlab,max(abs(lx-vx(:))),max(abs(ly-vy(:)))]);
            dmass=max(dmass,max(abs(lm-gather(massbin(:)))));
        end
    end
    fprintf('AutoCorr brute force B, assets TransitionProbs (type %s) PType vs binned dense kernel (all ages), this should be zero: %.3e, NaN-pattern mismatches %i \n',nm,dtp,nanmm)
    fprintf('AutoCorr brute force B, assets TransitionProbs (type %s) row sums minus one (bins with mass, all ages), this should be zero: %.3e \n',nm,drow)
    fprintf('AutoCorr brute force B, assets TransitionValues_j and _jplus1 (type %s) vs the unique values (all ages), this should be zero: %.3e \n',nm,dlab)
    fprintf('AutoCorr brute force B, assets TransitionMass_j (type %s) vs the origin bin masses (all ages), this should be zero: %.3e \n',nm,dmass)
end

% Grouped vs by-hand pooling of the two types' brute-force joint distributions
Dist_u_solo=struct();
for ii=1:N_i
    nm=Names_i{ii};
    Dist_u_solo.(nm)=gpuArray(reshape(gather(Dist_solo.(nm)),[N_s,N_j]));
end

% Grouped TransitionProbs (assets, earnings) vs a by-hand pooling of the types' dense-kernel joints. At age jj the pooled joint
% of (value at jj, value at jj+1) over the union of the types' value bins is sum over types of ptw(ii)*(type mass at jj) times
% the type's within-age joint (dist_ii(:,jj) normalised, times T_ii{jj}) binned by value; grouped row b is that divided by the
% pooled origin-bin mass (zero row if none), TransitionValues_j/_jplus1 are the union values and TransitionMass_j the pooled
% within-age shares. The grouped quintile TransitionProbs coarsen the same pooled joint to the quintile bins of the pooled
% age-jj value distribution (origin) and of the pooled age-(jj+1) value distribution (destination, weights
% ptw(ii)*(type mass at jj+1), from the types' actual age-(jj+1) distributions), bin q being the values up to the first sorted
% value whose cumulative mass exceeds q/5 (the single-type command's rule; part 20 checks that rule against the z model).
nq=5;
simoptions_PT_q=simoptions_PT;
simoptions_PT_q.transprobquantiles=nq;
ACP_PT_q=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_PT_q);
for fn={'assets','earnings'}
    ff=find(strcmp(FnNames,fn{1}));
    dtp=0; dlab=0; dmass=0; nbad=0; dq=0; nanq=0;
    for jj=1:N_j-1
        wj=zeros(N_i,1); wjp=zeros(N_i,1);
        for ii=1:N_i
            nm=Names_i{ii};
            wj(ii)=ptw(ii)*gather(sum(Dist_u_solo.(nm)(:,jj)));
            wjp(ii)=ptw(ii)*gather(sum(Dist_u_solo.(nm)(:,jj+1)));
        end
        if sum(wj)==0
            nbad=nbad+(~isempty(ACP_PT.(fn{1}).TransitionProbs{jj}));
            continue
        end
        xcell=cell(N_i,1); ycell=cell(N_i,1);
        for ii=1:N_i
            if wj(ii)>0
                xcell{ii}=gather(Vals_solo.(Names_i{ii}){ff}(:,jj));
                ycell{ii}=gather(Vals_solo.(Names_i{ii}){ff}(:,jj+1));
            end
        end
        vals_j=unique(vertcat(xcell{:})); vals_jp=unique(vertcat(ycell{:}));
        J_G=zeros(length(vals_j),length(vals_jp)); mass_j=zeros(length(vals_j),1); mass_jp=zeros(length(vals_jp),1);
        for ii=1:N_i
            if wj(ii)>0
                nm=Names_i{ii};
                [~,ix]=ismember(xcell{ii},vals_j); [~,iy]=ismember(ycell{ii},vals_jp);
                Sx=gpuArray(full(sparse(1:N_s,ix,1,N_s,length(vals_j)))); Sy=gpuArray(full(sparse(1:N_s,iy,1,N_s,length(vals_jp))));
                distj=Dist_u_solo.(nm)(:,jj)/sum(Dist_u_solo.(nm)(:,jj));
                distjp=Dist_u_solo.(nm)(:,jj+1)/sum(Dist_u_solo.(nm)(:,jj+1));
                J_G=J_G+wj(ii)*gather(Sx'*(distj.*T_solo.(nm){jj})*Sy);
                mass_j=mass_j+wj(ii)*gather(Sx'*distj);
                mass_jp=mass_jp+wjp(ii)*gather(Sy'*distjp);
            end
        end
        mass_safe=mass_j; mass_safe(mass_j==0)=1;
        P_bf=J_G./mass_safe;
        P=gather(ACP_PT.(fn{1}).TransitionProbs{jj});
        lx=gather(ACP_PT.(fn{1}).TransitionValues_j{jj}(:)); ly=gather(ACP_PT.(fn{1}).TransitionValues_jplus1{jj}(:)); lm=gather(ACP_PT.(fn{1}).TransitionMass_j{jj}(:));
        if ~isequal(size(P),size(P_bf)) || numel(lx)~=numel(vals_j) || numel(ly)~=numel(vals_jp) || numel(lm)~=numel(vals_j)
            nbad=nbad+1;
        else
            dtp=max(dtp,max(abs(P(:)-P_bf(:))));
            dlab=max([dlab,max(abs(lx-vals_j(:))),max(abs(ly-vals_jp(:)))]);
            dmass=max(dmass,max(abs(lm-mass_j(:)/sum(wj))));
        end
        if sum(wjp)>0
            cum_j=cumsum(mass_j)/sum(wj); idx_j=ones(length(vals_j),1);
            for qq=1:nq-1
                cutoff=vals_j(find(cum_j>qq/nq,1,'first')); idx_j(vals_j>cutoff)=qq+1;
            end
            cum_jp=cumsum(mass_jp)/sum(wjp); idx_jp=ones(length(vals_jp),1);
            for qq=1:nq-1
                cutoff=vals_jp(find(cum_jp>qq/nq,1,'first')); idx_jp(vals_jp>cutoff)=qq+1;
            end
            Sq_j=full(sparse(1:length(vals_j),idx_j,1,length(vals_j),nq)); Sq_jp=full(sparse(1:length(vals_jp),idx_jp,1,length(vals_jp),nq));
            massq=Sq_j'*mass_j; massq(massq==0)=1;
            Pq_bf=(Sq_j'*J_G*Sq_jp)./massq;
            Pq=gather(ACP_PT_q.(fn{1}).TransitionProbs(:,:,jj));
            dd=abs(Pq-Pq_bf); dd=dd(~isnan(dd)); dq=max(dq,max([dd(:);0])); nanq=nanq+sum(isnan(Pq(:))~=isnan(Pq_bf(:)));
        end
    end
    fprintf('AutoCorr brute force B, %s grouped TransitionProbs vs by-hand pooling of the dense-kernel joints (all ages), this should be zero: %.3e, wrongly sized cells %i \n',fn{1},dtp,nbad)
    fprintf('AutoCorr brute force B, %s grouped TransitionValues_j and _jplus1 vs the union of the types'' values (all ages), this should be zero: %.3e \n',fn{1},dlab)
    fprintf('AutoCorr brute force B, %s grouped TransitionMass_j vs the pooled origin bin masses (all ages), this should be zero: %.3e \n',fn{1},dmass)
    fprintf('AutoCorr brute force B, %s grouped quintile TransitionProbs vs by-hand pooling (all ages), this should be zero: %.3e, NaN-pattern mismatches %i \n',fn{1},dq,nanq)
    fprintf('AutoCorr brute force B, %s grouped quintile TransitionProbs shape [%i %i %i], expected [%i %i %i] \n',fn{1},size(ACP_PT_q.(fn{1}).TransitionProbs,1),size(ACP_PT_q.(fn{1}).TransitionProbs,2),size(ACP_PT_q.(fn{1}).TransitionProbs,3),nq,nq,N_j-1)
end
for ff=1:length(FnNames)
    fn=FnNames{ff};
    % Grouped Mean and StdDeviation at each age, by hand: weight ptw(ii)*(mass of type ii at age jj)
    MeanG_bf=nan(1,N_j); StdG_bf=nan(1,N_j);
    MeanRG_bf=nan(length(RestNames),N_j); StdRG_bf=nan(length(RestNames),N_j);
    for jj=1:N_j
        wmass=zeros(N_i,1); mx=zeros(N_i,1); mxx=zeros(N_i,1);
        wmassR=zeros(N_i,length(RestNames)); mxR=zeros(N_i,length(RestNames)); mxxR=zeros(N_i,length(RestNames));
        for ii=1:N_i
            nm=Names_i{ii};
            x=Vals_solo.(nm){ff}(:,jj);
            wmass(ii)=gather(ptw(ii)*sum(Dist_u_solo.(nm)(:,jj)));
            mx(ii)=gather(ptw(ii)*sum(Dist_u_solo.(nm)(:,jj).*x));
            mxx(ii)=gather(ptw(ii)*sum(Dist_u_solo.(nm)(:,jj).*x.^2));
            for rr=1:length(RestNames)
                mr=Dist_u_solo.(nm)(:,jj).*Rvals_solo.(nm){rr}(:,jj);
                wmassR(ii,rr)=gather(ptw(ii)*sum(mr));
                mxR(ii,rr)=gather(ptw(ii)*sum(mr.*x));
                mxxR(ii,rr)=gather(ptw(ii)*sum(mr.*x.^2));
            end
        end
        if sum(wmass)>0
            MeanG_bf(jj)=sum(mx)/sum(wmass);
            StdG_bf(jj)=sqrt(max(sum(mxx)/sum(wmass)-MeanG_bf(jj)^2,0));
        end
        for rr=1:length(RestNames)
            if sum(wmassR(:,rr))>0
                MeanRG_bf(rr,jj)=sum(mxR(:,rr))/sum(wmassR(:,rr));
                StdRG_bf(rr,jj)=sqrt(max(sum(mxxR(:,rr))/sum(wmassR(:,rr))-MeanRG_bf(rr,jj)^2,0));
            end
        end
    end
    x=gather(ACP_PT.(fn).Mean(:)); y=MeanG_bf(:); dd=abs(x-y); dd=dd(~isnan(dd));
    fprintf('AutoCorr brute force B, %s grouped Mean vs by-hand pooling, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,max([dd;0]),sum(isnan(x)~=isnan(y)))
    x=gather(ACP_PT.(fn).StdDeviation(:)); y=StdG_bf(:); dd=abs(x-y); dd=dd(~isnan(dd));
    fprintf('AutoCorr brute force B, %s grouped StdDeviation vs by-hand pooling, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,max([dd;0]),sum(isnan(x)~=isnan(y)))
    for rr=1:length(RestNames)
        rn=RestNames{rr};
        x=gather(ACP_PT.(rn).(fn).Mean(:)); y=MeanRG_bf(rr,:)'; dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('AutoCorr brute force B, %s %s grouped Mean vs by-hand pooling, this should be zero: %.3e, NaN-pattern mismatches %i \n',rn,fn,max([dd;0]),sum(isnan(x)~=isnan(y)))
        x=gather(ACP_PT.(rn).(fn).StdDeviation(:)); y=StdRG_bf(rr,:)'; dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('AutoCorr brute force B, %s %s grouped StdDeviation vs by-hand pooling, this should be zero: %.3e, NaN-pattern mismatches %i \n',rn,fn,max([dd;0]),sum(isnan(x)~=isnan(y)))
    end

    % Grouped pairs at each horizon, by hand: pool the types' unnormalised joint distributions with weight ptw(ii)
    for hh=1:length(horizons_PT)
        kk=horizons_PT(hh);
        AutoCovG_bf=nan(1,N_j-kk); AutoCorrG_bf=nan(1,N_j-kk);
        AutoCovRG_bf=nan(length(RestNames),N_j-kk); AutoCorrRG_bf=nan(length(RestNames),N_j-kk); PairMassG_bf=nan(length(RestNames),N_j-kk);
        PairMeanxG_bf=nan(length(RestNames),N_j-kk); PairMeanyG_bf=nan(length(RestNames),N_j-kk); PairStdxG_bf=nan(length(RestNames),N_j-kk); PairStdyG_bf=nan(length(RestNames),N_j-kk);
        for j0=1:N_j-kk
            Ju=cell(N_i,1); xs=cell(N_i,1); ys=cell(N_i,1);
            for ii=1:N_i
                nm=Names_i{ii};
                Pk=eye(N_s,'gpuArray');
                for s2=0:kk-1
                    Pk=Pk*T_solo.(nm){j0+s2};
                end
                Ju{ii}=ptw(ii)*(Dist_u_solo.(nm)(:,j0).*Pk); % unnormalised joint of type ii, weighted by its ptweight
                xs{ii}=Vals_solo.(nm){ff}(:,j0);
                ys{ii}=Vals_solo.(nm){ff}(:,j0+kk);
            end
            % unrestricted
            total=0; sx=0; sy=0;
            for ii=1:N_i
                total=total+sum(Ju{ii}(:)); sx=sx+sum(sum(Ju{ii},2).*xs{ii}); sy=sy+sum(sum(Ju{ii},1)'.*ys{ii});
            end
            if total>0
                mux=sx/total; muy=sy/total;
                vx=0; vy=0; cv=0;
                for ii=1:N_i
                    vx=vx+sum(sum(Ju{ii},2).*(xs{ii}-mux).^2); vy=vy+sum(sum(Ju{ii},1)'.*(ys{ii}-muy).^2);
                    cv=cv+sum(sum(Ju{ii}.*((xs{ii}-mux)*(ys{ii}-muy)')));
                end
                AutoCovG_bf(j0)=gather(cv/total);
                if sqrt(vx*vy)/total>1e-15
                    AutoCorrG_bf(j0)=gather(cv/sqrt(vx*vy));
                end
            end
            % restricted
            for rr=1:length(RestNames)
                total=0; sx=0; sy=0; Jr=cell(N_i,1);
                for ii=1:N_i
                    nm=Names_i{ii};
                    Jr{ii}=Ju{ii}.*(Rvals_solo.(nm){rr}(:,j0)*Rvals_solo.(nm){rr}(:,j0+kk)');
                    total=total+sum(Jr{ii}(:)); sx=sx+sum(sum(Jr{ii},2).*xs{ii}); sy=sy+sum(sum(Jr{ii},1)'.*ys{ii});
                end
                PairMassG_bf(rr,j0)=gather(total);
                if total>0
                    mux=sx/total; muy=sy/total;
                    vx=0; vy=0; cv=0;
                    for ii=1:N_i
                        vx=vx+sum(sum(Jr{ii},2).*(xs{ii}-mux).^2); vy=vy+sum(sum(Jr{ii},1)'.*(ys{ii}-muy).^2);
                        cv=cv+sum(sum(Jr{ii}.*((xs{ii}-mux)*(ys{ii}-muy)')));
                    end
                    PairMeanxG_bf(rr,j0)=gather(mux); PairMeanyG_bf(rr,j0)=gather(muy);
                    PairStdxG_bf(rr,j0)=gather(sqrt(vx/total)); PairStdyG_bf(rr,j0)=gather(sqrt(vy/total));
                    AutoCovRG_bf(rr,j0)=gather(cv/total);
                    if sqrt(vx*vy)/total>1e-15
                        AutoCorrRG_bf(rr,j0)=gather(cv/sqrt(vx*vy));
                    end
                end
            end
        end
        x=gather(ACP_PT.(fn).(['AutoCovariance',hstr_PT{hh}])(:)); y=AutoCovG_bf(:); dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('AutoCorr brute force B, %s grouped AutoCovariance%s vs by-hand pooling, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,hstr_PT{hh},max([dd;0]),sum(isnan(x)~=isnan(y)))
        x=gather(ACP_PT.(fn).(['AutoCorrelation',hstr_PT{hh}])(:)); y=AutoCorrG_bf(:); dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('AutoCorr brute force B, %s grouped AutoCorrelation%s vs by-hand pooling, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,hstr_PT{hh},max([dd;0]),sum(isnan(x)~=isnan(y)))
        for rr=1:length(RestNames)
            rn=RestNames{rr};
            bfnames={'AutoCovariance','AutoCorrelation','PairMass','PairMean_j','PairMean_jplusk','PairStdDeviation_j','PairStdDeviation_jplusk'};
            bfvals={AutoCovRG_bf,AutoCorrRG_bf,PairMassG_bf,PairMeanxG_bf,PairMeanyG_bf,PairStdxG_bf,PairStdyG_bf};
            for bb=1:length(bfnames)
                x=gather(ACP_PT.(rn).(fn).([bfnames{bb},hstr_PT{hh}])(:)); y=bfvals{bb}(rr,:)'; dd=abs(x-y); dd=dd(~isnan(dd));
                fprintf('AutoCorr brute force B, %s %s grouped %s%s vs by-hand pooling, this should be zero: %.3e, NaN-pattern mismatches %i \n',rn,fn,bfnames{bb},hstr_PT{hh},max([dd;0]),sum(isnan(x)~=isnan(y)))
            end
        end
    end
end
for rr=1:length(RestNames)
    rn=RestNames{rr};
    ByAge_bf=zeros(1,N_j);
    for ii=1:N_i
        nm=Names_i{ii};
        ByAge_bf=ByAge_bf+ptw(ii)*gather(sum(Dist_u_solo.(nm).*Rvals_solo.(nm){rr},1));
    end
    fprintf('AutoCorr brute force B, %s RestrictedSampleMass.ByAge vs by hand, this should be zero: %.3e \n',rn,max(abs(gather(ACP_PT.(rn).RestrictedSampleMass.ByAge(:))-ByAge_bf(:))))
    fprintf('AutoCorr brute force B, %s RestrictedSampleMass.Total vs by hand, this should be zero: %.3e \n',rn,abs(gather(ACP_PT.(rn).RestrictedSampleMass.Total)-sum(ByAge_bf)))
end

% A type of zero weight: grouped equals the other type alone
Params_PT0=Params_PT;
Params_PT0.(PTypeDistParamNames{1})=[0;1];
StationaryDist_PT0=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_PT0,simoptions);
ACP_PT0=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_PT0,Policy_PT,FnsToEvaluate_PT,Params_PT0,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_PT);
for ff=1:length(FnNames)
    fn=FnNames{ff};
    statnames={'Mean','StdDeviation'};
    for hh=1:length(horizons_PT)
        statnames=[statnames,{['AutoCovariance',hstr_PT{hh}],['AutoCorrelation',hstr_PT{hh}]}]; %#ok<AGROW>
    end
    for ss=1:length(statnames)
        x=gather(ACP_PT0.(fn).(statnames{ss})(:)); y=gather(ACP_solo.high.(fn).(statnames{ss})(:)); dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('AutoCorr brute force B, zero-weight type, %s grouped %s vs the other type alone, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,statnames{ss},max([dd;0]),sum(isnan(x)~=isnan(y)))
    end
    for rr=1:length(RestNames)
        rn=RestNames{rr};
        statnames=fieldnames(ACP_solo.high.(rn).(fn));
        for ss=1:length(statnames)
            x=gather(ACP_PT0.(rn).(fn).(statnames{ss})(:)); y=gather(ACP_solo.high.(rn).(fn).(statnames{ss})(:)); dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('AutoCorr brute force B, zero-weight type, %s %s grouped %s vs the other type alone, this should be zero: %.3e, NaN-pattern mismatches %i \n',rn,fn,statnames{ss},max([dd;0]),sum(isnan(x)~=isnan(y)))
        end
    end
end

output=struct();

end
