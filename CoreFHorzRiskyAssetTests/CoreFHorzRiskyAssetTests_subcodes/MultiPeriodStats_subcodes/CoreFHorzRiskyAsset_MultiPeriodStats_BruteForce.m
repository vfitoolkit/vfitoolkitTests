function output=CoreFHorzRiskyAsset_MultiPeriodStats_BruteForce(n_d1,n_d2,n_d3,n_d4,d1_grid,d2_grid,d3_grid,d4_grid,n_a,a_grid,n_a1,n_a2,a1_grid,a2_grid,n_z,z_grid,pi_z,N_j,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline)
% riskyasset: the multi-period statistics with a risky asset (a2'=aprimeFn(riskyshare,savings,u)=(1+r)*(1-riskyshare)*savings
% +(1+r+u)*riskyshare*savings, off the grid, so every state goes to two a2' grid points for each u, each with probability pi_u(u)
% times the interpolation weight; a2' does not depend on a2) against a brute force. The same as
% CoreFHorzExpAssetU_MultiPeriodStats_BruteForce in the ExpAssetU bank, for the risky asset. For each age jj the one-period
% transition matrix T_jj(s,s') is built here, from the Policy indexes and the law of motion of a2 written out again (not the
% toolkit's CreateaprimePolicyRiskyAsset): the linear-interpolation weights of a2' on the a2 grid (clipped at the ends of the grid),
% the second-layer weights of a1' with gridinterplayer=1, the semi-exogenous transition from the SemiExoStateFn itself, then pi_z
% and pi_e.
% Checks:
%   push consistency: StationaryDist(:,jj)'*T_jj equals StationaryDist(:,jj+1), which ties this T to the agent distribution;
%   EvalFnOnAgentDist_AutoCorrTransProbs_FHorz: Mean, StdDeviation, AutoCovariance/AutoCorrelation (horizon 1) and
%       AutoCovariance_horizon2 of each function, and AutoCovariance_horizon2 of its 2-period window, from T;
%   LifeCycleProfiles_FHorz_Case1: Mean and Variance of the 2-period window of each function, from T;
%   EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz: CovarianceMatrix_2periodAverage, from T.
% E[x_a x_b] (a<b) is mu_a'(x_a.*(T_a*...*T_{b-1}*x_b))/m_a, with mu_a the agent distribution at age a and m_a its mass.
% Models: (1) d1, a1 (safe asset), the risky asset, z and e (refine_d=[1,1,1]); (2) the same with gridinterplayer=1; (3) d1, a1,
% the risky asset, semiz, z and e (refine_d=[1,1,1,1], RiskyAssetSemiExo); (4) no d1, no a1, no z, no e (refine_d=[0,1,1],
% RiskyAssetu_noz).

fprintf('\n========== CoreFHorzRiskyAsset_MultiPeriodStats_BruteForce ==========\n')

pi_e_cpu=gather(vfoptionsbaseline.pi_e(:));
n_e=vfoptionsbaseline.n_e;
u_cpu=gather(vfoptionsbaseline.u_grid(:)); pi_u_cpu=gather(vfoptionsbaseline.pi_u(:)); N_u=length(u_cpu);

for mm=1:4
    vfoptions=struct();
    vfoptions.riskyasset=1;
    vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
    vfoptions.n_u=vfoptionsbaseline.n_u; vfoptions.u_grid=vfoptionsbaseline.u_grid; vfoptions.pi_u=vfoptionsbaseline.pi_u;
    if mm==1 || mm==2
        modelstr=['d1 a1 z e',repmat(' gridinterplayer',1,mm==2)];
        n_d=[n_d1,n_d2,n_d3]; d_grid=[d1_grid;d2_grid;d3_grid];
        vfoptions.refine_d=[1,1,1];
        n_am=[n_a1,n_a2]; a_gridm=[a1_grid;a2_grid];
        n_zm=n_z; z_gridm=z_grid; pi_zm=pi_z;
        vfoptions.n_e=n_e; vfoptions.e_grid=vfoptionsbaseline.e_grid; vfoptions.pi_e=vfoptionsbaseline.pi_e;
        if mm==2
            vfoptions.gridinterplayer=1;
            vfoptions.ngridinterp=5;
        end
        ReturnFn=@(h,savings,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,eta,varphi,r_a1,agej,Jr,pension) ReturnFn_d1_z_e_nosemiz_withA1(h,savings,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,eta,varphi,r_a1,agej,Jr,pension);
        FnsToEvaluate=struct();
        FnsToEvaluate.safeasset=@(h,riskyshare,savings,a1prime,a1,a2,z,e) a1;
        FnsToEvaluate.riskyasset=@(h,riskyshare,savings,a1prime,a1,a2,z,e) a2;
        FnsToEvaluate.earnings=@(h,riskyshare,savings,a1prime,a1,a2,z,e,w,kappa_j) w*kappa_j*h*z*e;
        jequaloneDist=zeros([n_am,n_zm,n_e],'gpuArray');
        jequaloneDist(1,ceil(n_a2/2),ceil(n_zm/2),ceil(n_e/2))=1;
        n_a1m=n_a1; n_a2m=n_a2;
        N_semizr=1; N_zr=n_zm; N_er=n_e;
        rowd2=2; rowd3=3; rowd4=0; rowa1p=4;
    elseif mm==3
        modelstr='d1 a1 semiz z e';
        n_d=[n_d1,n_d2,n_d3,n_d4]; d_grid=[d1_grid;d2_grid;d3_grid;d4_grid];
        vfoptions.refine_d=[1,1,1,1];
        n_am=[n_a1,n_a2]; a_gridm=[a1_grid;a2_grid];
        n_zm=n_z; z_gridm=z_grid; pi_zm=pi_z;
        vfoptions.n_e=n_e; vfoptions.e_grid=vfoptionsbaseline.e_grid; vfoptions.pi_e=vfoptionsbaseline.pi_e;
        vfoptions.n_semiz=vfoptionsbaseline.n_semiz; vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid; vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
        ReturnFn=@(h,savings,dsemiz,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,eta,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_d1_z_e_semiz_withA1(h,savings,dsemiz,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,eta,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
        FnsToEvaluate=struct();
        FnsToEvaluate.safeasset=@(h,riskyshare,savings,dsemiz,a1prime,a1,a2,semiz,z,e) a1;
        FnsToEvaluate.riskyasset=@(h,riskyshare,savings,dsemiz,a1prime,a1,a2,semiz,z,e) a2;
        FnsToEvaluate.earnings=@(h,riskyshare,savings,dsemiz,a1prime,a1,a2,semiz,z,e,w,kappa_j) w*kappa_j*h*semiz*z*e;
        jequaloneDist=zeros([n_am,vfoptions.n_semiz,n_zm,n_e],'gpuArray');
        jequaloneDist(1,ceil(n_a2/2),2,ceil(n_zm/2),ceil(n_e/2))=1;
        n_a1m=n_a1; n_a2m=n_a2;
        N_semizr=vfoptions.n_semiz; N_zr=n_zm; N_er=n_e;
        rowd2=2; rowd3=3; rowd4=4; rowa1p=5;
        semiz_grid=gather(vfoptions.semiz_grid(:));
    else
        modelstr='noa1 noz noe';
        n_d=[n_d2,n_d3]; d_grid=[d2_grid;d3_grid];
        vfoptions.refine_d=[0,1,1];
        n_am=n_a; a_gridm=a_grid;
        n_zm=0; z_gridm=[]; pi_zm=[];
        ReturnFn=@(savings,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz(savings,a,r,w,kappa_j,sigma,agej,Jr,pension);
        FnsToEvaluate=struct();
        FnsToEvaluate.riskyasset=@(riskyshare,savings,a) a;
        FnsToEvaluate.savings=@(riskyshare,savings,a) savings;
        FnsToEvaluate.riskysavings=@(riskyshare,savings,a) riskyshare*savings;
        jequaloneDist=zeros([n_am,1],'gpuArray');
        jequaloneDist(ceil(n_am/2))=1; % mid-grid wealth, so that the agents have savings to put in the risky asset from age 1
        n_a1m=1; n_a2m=n_a; % no standard endogenous state
        N_semizr=1; N_zr=1; N_er=1;
        rowd2=1; rowd3=2; rowd4=0; rowa1p=0;
    end
    a2g=gather(a_gridm(end-n_a2m+1:end)); a2g=a2g(:);
    if n_zm>0
        pi_z_cpu=gather(pi_zm);
    else
        pi_z_cpu=1;
    end
    simoptions=vfoptions;
    simoptions.d_grid=d_grid;
    simoptions.a_grid=a_gridm;
    if n_zm>0
        simoptions.pi_z=pi_zm; % the M>1 windows of the cross-section commands and LifeCycleProfiles need the z transitions
    end
    FnNames=fieldnames(FnsToEvaluate);
    nF=length(FnNames);

    [~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_am,n_zm,N_j,d_grid,a_gridm,z_gridm,pi_zm,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_am,n_zm,N_j,pi_zm,Params,simoptions);

    soA=simoptions; soA.timehorizons=[1,2]; soA.ntimeperiods=[1,2];
    AC=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_am,n_zm,N_j,d_grid,a_gridm,z_gridm,pi_zm,soA);
    soC=simoptions; soC.ntimeperiods=[1,2];
    LCP=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_am,n_zm,N_j,d_grid,a_gridm,z_gridm,soC);
    ACS=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_am,n_zm,N_j,d_grid,a_gridm,z_gridm,soC);
    ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy,FnsToEvaluate,Params,[],n_d,n_am,n_zm,N_j,d_grid,a_gridm,z_gridm,simoptions);

    %% The brute-force transition matrices
    N_a=n_a1m*n_a2m;
    N_s=N_a*N_semizr*N_zr*N_er;
    mu=gather(reshape(StationaryDist,[N_s,N_j]));
    X=cell(nF,1);
    for ff=1:nF
        X{ff}=gather(reshape(ValuesOnGrid.(FnNames{ff}),[N_s,N_j]));
    end
    PolicyK=gather(reshape(Policy,[size(Policy,1),N_s,N_j]));
    ss=(1:1:N_s)';
    rest=floor((ss-1)/N_a);
    isz=mod(rest,N_semizr)+1;
    iz=mod(floor(rest/N_semizr),N_zr)+1;
    T=cell(N_j-1,1);
    for jj=1:N_j-1
        rsval=d2_grid(PolicyK(rowd2,:,jj)'); rsval=rsval(:); % riskyshare
        svval=d3_grid(PolicyK(rowd3,:,jj)'); svval=svval(:); % savings
        % a2' and its corners, for each u with probability pi_u (written out again: the law of motion of the risky asset and the
        % linear interpolation, with everything off the bottom of the grid on the bottom point and off the top on the top point)
        corneridx=zeros(N_s,0); cornerprob=zeros(N_s,0);
        for uu=1:N_u
            a2p=(1+Params.r)*svval; % negative or zero savings are all safe
            pos=(svval>0);
            a2p(pos)=(1+Params.r)*(1-rsval(pos)).*svval(pos)+(1+Params.r+u_cpu(uu))*rsval(pos).*svval(pos);
            lo=min(max(sum(a2g'<=a2p,2),1),n_a2m-1);
            wl=min(max((a2g(lo+1)-a2p)./(a2g(lo+1)-a2g(lo)),0),1);
            corneridx=[corneridx,lo,lo+1]; %#ok<AGROW> % the a2 index
            cornerprob=[cornerprob,wl*pi_u_cpu(uu),(1-wl)*pi_u_cpu(uu)]; %#ok<AGROW>
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
        % semiz' (from the SemiExoStateFn, at the d4 of the state)
        if N_semizr>1
            d4val=d4_grid(PolicyK(rowd4,:,jj)'); d4val=d4val(:);
            psz=zeros(N_s,N_semizr);
            for s2=1:N_semizr
                psz(:,s2)=arrayfun(@(nn,dd) CoreFHorzRiskyAssetSetup_SemiExoStateFn(nn,semiz_grid(s2),dd,Params.probfindjob,Params.problosejob), semiz_grid(isz), d4val);
            end
        else
            psz=ones(N_s,1);
        end
        rows=zeros(0,1); cols=zeros(0,1); vals=zeros(0,1);
        for c1=1:size(a1pidx,2)
            for c2=1:size(corneridx,2)
                for s2=1:N_semizr
                    for z2=1:N_zr
                        for e2=1:N_er
                            dest=a1pidx(:,c1)+n_a1m*(corneridx(:,c2)-1)+N_a*(s2-1)+N_a*N_semizr*(z2-1)+N_a*N_semizr*N_zr*(e2-1);
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
    fprintf('MultiPeriodStats riskyasset brute force (%s): push consistency, StationaryDist(:,jj)''*T_jj vs StationaryDist(:,jj+1), this should be zero: %.3e \n',modelstr,pushdiff)

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
        fprintf('MultiPeriodStats riskyasset brute force (%s): %s vs brute force, this should be zero: %.3e \n',modelstr,labels{ll},maxd(ll))
        fprintf('MultiPeriodStats riskyasset brute force (%s): %s, entries NaN in only one (or of a different size), this should be zero: %i \n',modelstr,labels{ll},nnanbad(ll))
    end
    % Not vacuous: the risky asset really is off the grid (some state puts weight on both of its a2 corners), and u moves it (a
    % state goes to more than the two corners of one a2')
    fprintf('MultiPeriodStats riskyasset brute force (%s): largest weight on the less likely a2 corner, first u (should be well above zero): %.3f \n',modelstr,full(max(min(cornerprob(:,1),cornerprob(:,2))/pi_u_cpu(1))))
    fprintf('MultiPeriodStats riskyasset brute force (%s): number of distinct a2 destinations from one state across u (should be well above 2): %i \n',modelstr,max(arrayfun(@(s) numel(unique(corneridx(s,:))),1:N_s)))
    clear T Exy
end

output=struct();

end
