function output=CoreFHorzExpAssete_MultiPeriodStats_BruteForce(n_d1,n_d2,n_d3,d1_grid,d2_grid,d3_grid,n_a,a1_grid,a2_grid,n_a2_1,n_a2_2,a2_1_grid,a2_2_grid,n_z,z_grid,pi_z,N_j,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline)
% experienceassete: the multi-period statistics with a genuine experience asset that depends on e (a2'=phi1*(1-d2)*e+(1-phi2)*a2,
% with the e of the state, off the grid, so every state goes to two a2' grid points, or four with two experience assets) against
% a brute force. The same as CoreFHorzExpAsset_MultiPeriodStats_BruteForce in the ExpAsset bank (and _z in the ExpAssetz bank), with
% e in the law of motion, so a mistake in which e each state's a2' is computed with would show. For each age jj the
% one-period transition matrix T_jj(s,s') is built here, from the Policy indexes and the law of motion of a2 written out again (not
% the toolkit's CreateaprimePolicyExperienceAssete): the linear-interpolation weights of a2' on the a2 grid, the second-layer weights
% of a1' with gridinterplayer=1, the semi-exogenous transition from the SemiExoStateFn itself, then pi_z and pi_e.
% Checks:
%   push consistency: StationaryDist(:,jj)'*T_jj equals StationaryDist(:,jj+1), which ties this T to the agent distribution;
%   EvalFnOnAgentDist_AutoCorrTransProbs_FHorz: Mean, StdDeviation, AutoCovariance/AutoCorrelation (horizon 1) and
%       AutoCovariance_horizon2 of each function, and AutoCovariance_horizon2 of its 2-period window, from T;
%   LifeCycleProfiles_FHorz_Case1: Mean and Variance of the 2-period window of each function, from T;
%   EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz: CovarianceMatrix_2periodAverage, from T.
% E[x_a x_b] (a<b) is mu_a'(x_a.*(T_a*...*T_{b-1}*x_b))/m_a, with mu_a the agent distribution at age a and m_a its mass.
% The models are listed in the code (each has the shocks the experience asset needs). The a1 models use a coarser a1 grid (31 points,
% the same shape of grid): with semiz, z and e a state has up to 2*2*5*3 edges (twice that with gridinterplayer) and the brute force
% keeps one sparse T per age, which on the bank's 101-point a1 grid runs to GB. The check is the same on whatever grid is solved.

fprintf('\n========== CoreFHorzExpAssete_MultiPeriodStats_BruteForce ==========\n')

pi_z_cpu=gather(pi_z);
z_grid_cpu=gather(z_grid(:));
pi_e_cpu=gather(vfoptionsbaseline.pi_e(:));
n_e=vfoptionsbaseline.n_e;
e_grid_cpu=gather(vfoptionsbaseline.e_grid(:));
n_a1b=31;
a1_grid_b=5*linspace(0,1,n_a1b)'.^3;

for mm=1:4
    vfoptions=struct();
    vfoptions.n_e=n_e; vfoptions.e_grid=vfoptionsbaseline.e_grid; vfoptions.pi_e=vfoptionsbaseline.pi_e; % the experience asset depends on e, so every model has e
    if mm==1 || mm==2
        modelstr=['d1 a1 z e',repmat(' gridinterplayer',1,mm==2)];
        n_d=[n_d1,n_d2]; d_grid=[d1_grid;d2_grid];
        n_am=[n_a1b,n_a(2)]; a_grid=[a1_grid_b;a2_grid]; % the coarser a1 grid (see above)
        vfoptions.experienceassete=1;
        vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
        if mm==2
            vfoptions.gridinterplayer=1;
            vfoptions.ngridinterp=5;
        end
        ReturnFn=@(d1,d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_ExpAssete_d1_z_e_nosemiz(d1,d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(d1,d2,a1prime,a1,a2,z,e) a1;
        FnsToEvaluate.humancapital=@(d1,d2,a1prime,a1,a2,z,e) a2;
        FnsToEvaluate.earnings=@(d1,d2,a1prime,a1,a2,z,e,w,kappa_j) w*kappa_j*d1*d2*a2*z*e;
        jequaloneDist=zeros([n_am,n_z,n_e],'gpuArray');
        jequaloneDist(1,1,ceil(n_z/2),ceil(n_e/2))=1;
        n_a1=n_a1b; a1g=a1_grid_b; % the standard endogenous state
        n_a2v=n_a(2); % the experience asset(s)
        N_semizr=1; N_er=n_e;
        rowd2=2; rowd3=0; rowa1p=3;
    elseif mm==3
        modelstr='d1 a1 semiz z e';
        n_d=[n_d1,n_d2,n_d3]; d_grid=[d1_grid;d2_grid;d3_grid];
        n_am=[n_a1b,n_a(2)]; a_grid=[a1_grid_b;a2_grid]; % the coarser a1 grid (see above)
        vfoptions.n_semiz=vfoptionsbaseline.n_semiz; vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid; vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
        vfoptions.experienceassete=1;
        vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
        ReturnFn=@(d1,d2,d3,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssete_d1_z_e_semiz(d1,d2,d3,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(d1,d2,d3,a1prime,a1,a2,semiz,z,e) a1;
        FnsToEvaluate.humancapital=@(d1,d2,d3,a1prime,a1,a2,semiz,z,e) a2;
        FnsToEvaluate.earnings=@(d1,d2,d3,a1prime,a1,a2,semiz,z,e,w,kappa_j) w*kappa_j*d1*d2*a2*semiz*z*e;
        jequaloneDist=zeros([n_am,vfoptions.n_semiz,n_z,n_e],'gpuArray');
        jequaloneDist(1,2,1,ceil(n_z/2),ceil(n_e/2))=1;
        n_a1=n_a1b; a1g=a1_grid_b;
        n_a2v=n_a(2);
        N_semizr=vfoptions.n_semiz; N_er=n_e;
        rowd2=2; rowd3=3; rowa1p=4;
        semiz_grid=gather(vfoptions.semiz_grid(:));
    else
        modelstr='noa1 with2A2 z e';
        n_d=n_d2; d_grid=d2_grid;
        n_am=[n_a2_1,n_a2_2]; a_grid=[a2_1_grid;a2_2_grid];
        vfoptions.experienceassete=2;
        vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn_2A2;
        ReturnFn=@(d2,a2_1,a2_2,z,e,r,w,kappa_j,sigma,agej,Jr,pension,pensionrate) ReturnFn_ExpAssete_nod1_z_e_noa1_with2A2_nosemiz(d2,a2_1,a2_2,z,e,r,w,kappa_j,sigma,agej,Jr,pension,pensionrate);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(d2,a2_1,a2_2,z,e) a2_1;
        FnsToEvaluate.lifetimeearnings=@(d2,a2_1,a2_2,z,e) a2_2;
        FnsToEvaluate.earnings=@(d2,a2_1,a2_2,z,e,w,kappa_j) w*kappa_j*d2*a2_1*z*e;
        jequaloneDist=zeros([n_am,n_z,n_e],'gpuArray');
        jequaloneDist(1,1,ceil(n_z/2),ceil(n_e/2))=1;
        n_a1=1; a1g=[]; % no standard endogenous state
        n_a2v=[n_a2_1,n_a2_2];
        N_semizr=1; N_er=n_e;
        rowd2=1; rowd3=0; rowa1p=0;
    end
    simoptions=vfoptions;
    simoptions.d_grid=d_grid;
    simoptions.a_grid=a_grid;
    simoptions.pi_z=pi_z; % the M>1 windows of the cross-section commands and LifeCycleProfiles need the z transitions
    FnNames=fieldnames(FnsToEvaluate);
    nF=length(FnNames);

    [~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_am,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_am,n_z,N_j,pi_z,Params,simoptions);

    soA=simoptions; soA.timehorizons=[1,2]; soA.ntimeperiods=[1,2];
    AC=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_am,n_z,N_j,d_grid,a_grid,z_grid,pi_z,soA);
    soC=simoptions; soC.ntimeperiods=[1,2];
    LCP=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_am,n_z,N_j,d_grid,a_grid,z_grid,soC);
    ACS=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_am,n_z,N_j,d_grid,a_grid,z_grid,soC);
    ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy,FnsToEvaluate,Params,[],n_d,n_am,n_z,N_j,d_grid,a_grid,z_grid,simoptions);

    %% The brute-force transition matrices
    N_a=n_a1*prod(n_a2v);
    N_s=N_a*N_semizr*n_z*N_er;
    mu=gather(reshape(StationaryDist,[N_s,N_j]));
    X=cell(nF,1);
    for ff=1:nF
        X{ff}=gather(reshape(ValuesOnGrid.(FnNames{ff}),[N_s,N_j]));
    end
    PolicyK=gather(reshape(Policy,[size(Policy,1),N_s,N_j]));
    ss=(1:1:N_s)';
    ia=mod(ss-1,N_a)+1;
    ia1=mod(ia-1,n_a1)+1;
    ia2=floor((ia-1)/n_a1)+1; % the a2 index (linear over the a2 dims)
    rest=floor((ss-1)/N_a);
    isz=mod(rest,N_semizr)+1;
    iz=mod(floor(rest/N_semizr),n_z)+1;
    ie=floor(rest/(N_semizr*n_z))+1; % the e index (1 when there is no e)
    T=cell(N_j-1,1);
    for jj=1:N_j-1
        d2val=d2_grid(PolicyK(rowd2,:,jj)'); d2val=d2val(:);
        shockval=e_grid_cpu(ie); shockval=shockval(:); % the e of each state, which a2' depends on
        % a2' and its corners (written out again: the law of motion and the linear interpolation)
        if mm<=3
            a2val=a2_grid(ia2); a2val=a2val(:);
            a2p=Params.phi1*(1-d2val).*shockval+(1-Params.phi2)*a2val; % a2' depends on the e of the state
            lo=min(max(sum(a2_grid(:)'<=a2p,2),1),n_a2v-1);
            wl=min(max((a2_grid(lo+1)-a2p)./(a2_grid(lo+1)-a2_grid(lo)),0),1);
            corneridx=[lo,lo+1]; % [N_s,2], the a2 index
            cornerprob=[wl,1-wl];
        else
            ia2_1=mod(ia2-1,n_a2_1)+1;
            ia2_2=floor((ia2-1)/n_a2_1)+1;
            a21=a2_1_grid(ia2_1); a21=a21(:);
            a22=a2_2_grid(ia2_2); a22=a22(:);
            a21p=Params.phi1*(1-d2val).*shockval+(1-Params.phi2)*a21; % a2_1' depends on the e of the state
            a22p=Params.phi3*d2val.*a21+(1-Params.phi4)*a22;
            lo1=min(max(sum(a2_1_grid(:)'<=a21p,2),1),n_a2_1-1);
            wl1=min(max((a2_1_grid(lo1+1)-a21p)./(a2_1_grid(lo1+1)-a2_1_grid(lo1)),0),1);
            lo2=min(max(sum(a2_2_grid(:)'<=a22p,2),1),n_a2_2-1);
            wl2=min(max((a2_2_grid(lo2+1)-a22p)./(a2_2_grid(lo2+1)-a2_2_grid(lo2)),0),1);
            corneridx=[lo1+n_a2_1*(lo2-1), lo1+1+n_a2_1*(lo2-1), lo1+n_a2_1*lo2, lo1+1+n_a2_1*lo2];
            cornerprob=[wl1.*wl2, (1-wl1).*wl2, wl1.*(1-wl2), (1-wl1).*(1-wl2)];
        end
        % a1' (with the second-layer weights when gridinterplayer=1)
        if rowa1p==0
            a1pidx=ones(N_s,1); a1pprob=ones(N_s,1);
        elseif mm==2
            L2=PolicyK(end-1,:,jj)'; L2flag=PolicyK(end,:,jj)';
            L2(L2flag==1)=1;
            L2(L2flag==3)=vfoptions.ngridinterp+2;
            pu=(L2-1)/(vfoptions.ngridinterp+1);
            a1lo=PolicyK(rowa1p,:,jj)';
            a1pidx=[a1lo,a1lo+1]; a1pprob=[1-pu,pu];
        else
            a1pidx=PolicyK(rowa1p,:,jj)'; a1pprob=ones(N_s,1);
        end
        % semiz' (from the SemiExoStateFn, at the d3 of the state)
        if N_semizr>1
            d3val=d3_grid(PolicyK(rowd3,:,jj)'); d3val=d3val(:);
            psz=zeros(N_s,N_semizr);
            for s2=1:N_semizr
                psz(:,s2)=arrayfun(@(nn,dd) CoreFHorzExpAsseteSetup_SemiExoStateFn(nn,semiz_grid(s2),dd,Params.probfindjob,Params.problosejob), semiz_grid(isz), d3val);
            end
        else
            psz=ones(N_s,1);
        end
        rows=zeros(0,1); cols=zeros(0,1); vals=zeros(0,1);
        for c1=1:size(a1pidx,2)
            for c2=1:size(corneridx,2)
                for s2=1:N_semizr
                    for z2=1:n_z
                        for e2=1:N_er
                            dest=a1pidx(:,c1)+n_a1*(corneridx(:,c2)-1)+N_a*(s2-1)+N_a*N_semizr*(z2-1)+N_a*N_semizr*n_z*(e2-1);
                            pr=a1pprob(:,c1).*cornerprob(:,c2).*psz(:,s2).*pi_z_cpu(iz,z2);
                            if N_er>1
                                pr=pr*pi_e_cpu(e2);
                            end
                            rows=[rows;ss]; cols=[cols;dest]; vals=[vals;pr]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
        T{jj}=sparse(rows,cols,vals,N_s,N_s);
    end

    %% Push consistency
    pushdiff=0;
    for jj=1:N_j-1
        pushdiff=max(pushdiff,max(abs((mu(:,jj)'*T{jj})'-mu(:,jj+1))));
    end
    fprintf('MultiPeriodStats e brute force (%s): push consistency, StationaryDist(:,jj)''*T_jj vs StationaryDist(:,jj+1), this should be zero: %.3e \n',modelstr,pushdiff)

    %% The moments from T
    m=sum(mu,1); % the mass at each age
    Mean_bf=zeros(nF,N_j); Ex2_bf=zeros(nF,N_j); SD_bf=zeros(nF,N_j);
    for ff=1:nF
        Mean_bf(ff,:)=sum(mu.*X{ff},1)./m;
        Ex2_bf(ff,:)=sum(mu.*X{ff}.^2,1)./m;
        SD_bf(ff,:)=sqrt(sum(mu.*(X{ff}-Mean_bf(ff,:)).^2,1)./m); % centered (at age 1 all the mass is on one state, so this is exactly zero there)
    end
    % Exy{ff,gg}(a,b)=E[x_ff at age a times x_gg at age b], for b=a..a+3
    Exy=cell(nF,nF);
    for ff=1:nF
        for gg=1:nF
            Exy{ff,gg}=nan(N_j,N_j);
            for aa=1:N_j
                for bb=aa:min(aa+3,N_j)
                    v=X{gg}(:,bb);
                    for tt=bb-1:-1:aa
                        v=T{tt}*v;
                    end
                    Exy{ff,gg}(aa,bb)=sum(mu(:,aa).*X{ff}(:,aa).*v)/m(aa);
                end
            end
        end
    end

    maxd=zeros(1,7); nnanbad=zeros(1,7);
    labels={'AutoCorrTransProbs Mean','AutoCorrTransProbs StdDeviation','AutoCorrTransProbs AutoCovariance (horizon 1)','AutoCorrTransProbs AutoCorrelation (horizon 1)','AutoCorrTransProbs AutoCovariance_horizon2 and that of the 2-period window','LifeCycleProfiles 2-period window Mean and Variance','AgeConditionalStats_CrossSectionCovarCorr CovarianceMatrix_2periodAverage'};
    for ff=1:nF
        fn=FnNames{ff};
        cov1=zeros(1,N_j-1); cov2=zeros(1,N_j-2);
        for jj=1:N_j-1
            cov1(jj)=Exy{ff,ff}(jj,jj+1)-Mean_bf(ff,jj)*Mean_bf(ff,jj+1);
        end
        for jj=1:N_j-2
            cov2(jj)=Exy{ff,ff}(jj,jj+2)-Mean_bf(ff,jj)*Mean_bf(ff,jj+2);
        end
        denom=SD_bf(ff,1:end-1).*SD_bf(ff,2:end);
        corr1=cov1./denom;
        corr1(denom<=1e-15)=NaN; % as the command: no correlation where a standard deviation is zero
        % the 2-period window W_j=(x_j+x_{j+1})/2
        WMean=nan(1,N_j); WVar=nan(1,N_j); Wcov2=nan(1,N_j-2);
        for jj=1:N_j-1
            WMean(jj)=(Mean_bf(ff,jj)+Mean_bf(ff,jj+1))/2;
            WVar(jj)=(Ex2_bf(ff,jj)+2*Exy{ff,ff}(jj,jj+1)+Ex2_bf(ff,jj+1))/4-WMean(jj)^2;
        end
        for jj=1:N_j-3
            Wcov2(jj)=(Exy{ff,ff}(jj,jj+2)+Exy{ff,ff}(jj,jj+3)+Exy{ff,ff}(jj+1,jj+2)+Exy{ff,ff}(jj+1,jj+3))/4-WMean(jj)*WMean(jj+2);
        end
        Wname=[fn,'_2periodAverage'];
        cmp={AC.(fn).Mean,Mean_bf(ff,:),1; AC.(fn).StdDeviation,SD_bf(ff,:),2; AC.(fn).AutoCovariance,cov1,3; AC.(fn).AutoCorrelation,corr1,4; ...
            AC.(fn).AutoCovariance_horizon2,cov2,5; AC.(Wname).AutoCovariance_horizon2,Wcov2,5; ...
            LCP.(Wname).Mean,WMean,6; LCP.(Wname).Variance,WVar,6};
        for kk=1:size(cmp,1)
            xt=gather(double(cmp{kk,1}(:))); xb=cmp{kk,2}(:); ll=cmp{kk,3};
            if numel(xt)~=numel(xb)
                nnanbad(ll)=nnanbad(ll)+abs(numel(xt)-numel(xb))+1;
                continue
            end
            nnanbad(ll)=nnanbad(ll)+sum(isnan(xt)~=isnan(xb));
            both=~isnan(xt) & ~isnan(xb);
            if any(both)
                maxd(ll)=max(maxd(ll),max(abs(xt(both)-xb(both))));
            end
        end
        % the window covariances with every other function
        for gg=1:nF
            Wcov=nan(1,N_j);
            for jj=1:N_j-1
                EWW=(Exy{ff,gg}(jj,jj)+Exy{ff,gg}(jj,jj+1)+Exy{gg,ff}(jj,jj+1)+Exy{ff,gg}(jj+1,jj+1))/4;
                Wcov(jj)=EWW-(Mean_bf(ff,jj)+Mean_bf(ff,jj+1))/2*(Mean_bf(gg,jj)+Mean_bf(gg,jj+1))/2;
            end
            xt=gather(double(reshape(ACS.CovarianceMatrix_2periodAverage(ff,gg,:),[N_j,1]))); xb=Wcov(:);
            nnanbad(7)=nnanbad(7)+sum(isnan(xt)~=isnan(xb));
            both=~isnan(xt) & ~isnan(xb);
            if any(both)
                maxd(7)=max(maxd(7),max(abs(xt(both)-xb(both))));
            end
        end
    end
    for ll=1:7
        fprintf('MultiPeriodStats e brute force (%s): %s vs brute force, this should be zero: %.3e \n',modelstr,labels{ll},maxd(ll))
        fprintf('MultiPeriodStats e brute force (%s): %s, entries NaN in only one (or of a different size), this should be zero: %i \n',modelstr,labels{ll},nnanbad(ll))
    end
    % Not vacuous: the experience asset really is off the grid (some state puts weight on both of its a2 corners), and the
    % autocovariance of the experience asset is not zero
    fprintf('MultiPeriodStats e brute force (%s): largest weight on the less likely a2 corner (should be well above zero): %.3f \n',modelstr,full(max(min(cornerprob(:,1),1-cornerprob(:,1)))))
    clear T Exy
end

output=struct();

end
