function output=CoreInfHorz_AutoCorrBruteForce(n_a,n_z,n_e,a_grid,z_grid,pi_z,e_grid,pi_e,Params,DiscountFactorParamNames,n_d2_semiz,d2_grid_semiz,vfoptionsbaseline,simoptionsbaseline)
% EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz against a brute-force calculation written here independently of the
% command: Mean, StdDeviation, AutoCovariance and AutoCorrelation at horizons 1 and simoptions.timehorizons=[2,3,5]
% (the tperiodsK fields), and the TransitionProbs between unique values and between quintiles (simoptions.transprobquantiles=5)
% at the same horizons. (The FHorz command has the same brute force in CoreFHorzTests, TestFnsToEvaluate config 4.)
%
% Brute force: the one-period transition kernel T (N_states x N_states, dense) is built by loops from the Policy indexes and
% pi_z/pi_e (with the two-point grid-interpolation weights when gridinterplayer=1), the k-period kernel is T^k, and the joint
% distribution of (state now, state k periods on) is dist(s)*Tk(s,s'). Every moment is a plain sum over that joint
% distribution, and the transition probabilities between value bins are the joint distribution summed over the bins
% and divided by the origin bin's mass (NaN where the origin bin has no mass, which is what the command returns).
%
% Five models, each solved with and without the grid interpolation layer, so the command is exercised on
% z-only, e-only, z-and-e, semiz-and-z and semiz-z-and-e states (the e shock is an iid passed through simoptions.n_e/e_grid/pi_e, and the
% command keeps it apart from the markov z throughout).
%   z-only:  no d, markov z          (state (a,z),   ReturnFn_nod_z_noe_nosemiz)
%   e-only:  no d, iid e             (state (a,e),   ReturnFn_nod_noz_e_nosemiz)
%   z-and-e: no d, markov z, iid e   (state (a,z,e), ReturnFn_nod_z_e_nosemiz; z varies first, then e)
%   semiz-and-z: d2 only, semiz, markov z (state (a,semiz,z), ReturnFn_nod1_z_noe_semiz; semiz varies first, then z),
%            with and without the grid interpolation layer. The kernel uses pi_semiz(semiz,semiz'|d2) at the
%            policy d2, with pi_semiz built here from the SemiExoStateFn, and the shock is z (its transitions are pi_z^k
%            whatever the policy).
%   semiz-z-and-e: as semiz-and-z, with an iid e as well (state (a,semiz,z,e), ReturnFn_nod1_z_e_semiz; semiz varies
%            first, then z, then e), with and without the grid interpolation layer. The kernel also has pi_e(e'), and the shock is z*e (its transitions are
%            kron(pi_e,pi_z^k) whatever the policy, as in z-and-e).
% Checks, all exact identities (CoreSummary reads the diary back and applies its bar):
%   B1. the dense kernel's rows sum to one, and the stationary dist is a fixed point of it (close to zero: the dist is
%       only converged to the iteration tolerance)
%   B2. Mean, StdDeviation vs direct sums; AutoCovariance, AutoCorrelation at every horizon vs the brute force
%   B3. TransitionProbs between unique values (assets, and the shock itself) at every horizon vs the binned brute force,
%       rows summing to one where the origin bin has mass, and the shock's own TransitionProbs equal to pi_z^k (z) or to
%       pi_e' in every row (e)
%   B4. TransitionProbs between quintiles (assets) at every horizon vs the brute force with the same quantile definition
%       (bin q is the values up to the first sorted value whose cumulative mass exceeds q/5)
% Zero checks are printed as %.3e; NaN is compared as a pattern (count of entries NaN in one but not the other).

fprintf('\n========== CoreInfHorz_AutoCorrBruteForce ==========\n')

timehorizons=[2,3,5];
horizons=[1,timehorizons];
Kmax=max(horizons);
nq=5; % quintiles

pi_z_cpu=gather(pi_z);
pi_e_cpu=gather(pi_e(:));

for mm=1:5
    n_d=0; d_grid=[]; n_semizm=0;
    if mm==1
        mstr='z-only';
        n_zm=n_z; z_gridm=z_grid; pi_zm=pi_z; n_em=0;
        ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(aprime,a,z) a;
        FnsToEvaluate.earnings=@(aprime,a,z,w) w*z;
        FnsToEvaluate.consumption=@(aprime,a,z,r,w) (1+r)*a+w*z-aprime;
        FnsToEvaluate.shock=@(aprime,a,z) z;
    elseif mm==2
        mstr='e-only';
        n_zm=0; z_gridm=[]; pi_zm=[]; n_em=n_e;
        ReturnFn=@(aprime,a,e,r,w,sigma) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,sigma);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(aprime,a,e) a;
        FnsToEvaluate.earnings=@(aprime,a,e,w) w*e;
        FnsToEvaluate.consumption=@(aprime,a,e,r,w) (1+r)*a+w*e-aprime;
        FnsToEvaluate.shock=@(aprime,a,e) e;
    elseif mm==3
        mstr='z-and-e';
        n_zm=n_z; z_gridm=z_grid; pi_zm=pi_z; n_em=n_e;
        ReturnFn=@(aprime,a,z,e,r,w,sigma) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,sigma);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(aprime,a,z,e) a;
        FnsToEvaluate.earnings=@(aprime,a,z,e,w) w*z*e;
        FnsToEvaluate.consumption=@(aprime,a,z,e,r,w) (1+r)*a+w*z*e-aprime;
        FnsToEvaluate.shock=@(aprime,a,z,e) z*e; % the product takes n_z*n_e distinct values (the grids are exp of AR(1) grids, so no two products coincide)
    elseif mm==4
        mstr='semiz-and-z';
        n_zm=n_z; z_gridm=z_grid; pi_zm=pi_z; n_em=0;
        n_d=n_d2_semiz; d_grid=d2_grid_semiz; n_semizm=vfoptionsbaseline.n_semiz;
        ReturnFn=@(d2,aprime,a,semiz,z,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,sigma,uempbenefit,searcheffortcost);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(d2,aprime,a,semiz,z) a;
        FnsToEvaluate.earnings=@(d2,aprime,a,semiz,z,w) w*z*semiz;
        FnsToEvaluate.consumption=@(d2,aprime,a,semiz,z,r,w,uempbenefit) (1+r)*a+w*z*semiz+uempbenefit*(1-semiz)-aprime;
        FnsToEvaluate.shock=@(d2,aprime,a,semiz,z) z;
    else
        mstr='semiz-z-and-e';
        n_zm=n_z; z_gridm=z_grid; pi_zm=pi_z; n_em=n_e;
        n_d=n_d2_semiz; d_grid=d2_grid_semiz; n_semizm=vfoptionsbaseline.n_semiz;
        ReturnFn=@(d2,aprime,a,semiz,z,e,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_z_e_semiz(d2,aprime,a,semiz,z,e,r,w,sigma,uempbenefit,searcheffortcost);
        FnsToEvaluate=struct();
        FnsToEvaluate.assets=@(d2,aprime,a,semiz,z,e) a;
        FnsToEvaluate.earnings=@(d2,aprime,a,semiz,z,e,w) w*z*semiz*e;
        FnsToEvaluate.consumption=@(d2,aprime,a,semiz,z,e,r,w,uempbenefit) (1+r)*a+w*z*semiz*e+uempbenefit*(1-semiz)-aprime;
        FnsToEvaluate.shock=@(d2,aprime,a,semiz,z,e) z*e; % n_z*n_e distinct values, as in z-and-e
    end
    if n_semizm>0
        % pi_semiz(semiz,semiz',d2), straight from the SemiExoStateFn (on the grid values)
        semiz_gridm=vfoptionsbaseline.semiz_grid;
        pi_semiz_bf=zeros(n_semizm,n_semizm,n_d);
        for s1=1:n_semizm
            for s2=1:n_semizm
                for d2=1:n_d
                    pi_semiz_bf(s1,s2,d2)=CoreInfHorzSetup_SemiExoStateFn(semiz_gridm(s1),semiz_gridm(s2),d2_grid_semiz(d2),Params.probfindjob,Params.problosejob);
                end
            end
        end
    end
    FnNames=fieldnames(FnsToEvaluate);
    n_zr=max(n_zm,1); n_er=max(n_em,1); n_semizr=max(n_semizm,1);
    N_s=n_a*n_semizr*n_zr*n_er;
    vfoptions=struct();
    simoptions=struct();
    if n_em>0
        vfoptions.n_e=n_em; vfoptions.e_grid=e_grid; vfoptions.pi_e=pi_e;
        simoptions.n_e=n_em; simoptions.e_grid=e_grid; simoptions.pi_e=pi_e;
    end
    if n_semizm>0
        vfoptions.n_semiz=vfoptionsbaseline.n_semiz; vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid; vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
        simoptions.n_semiz=simoptionsbaseline.n_semiz; simoptions.semiz_grid=simoptionsbaseline.semiz_grid; simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
        simoptions.d_grid=d_grid;
    end

    for cfg=1:2
        if cfg==1
            cfgstr=[mstr,', no GI'];
            vfoptions_c=vfoptions; simoptions_c=simoptions; ngridinterp=0;
        else
            cfgstr=[mstr,', GI'];
            vfoptions_c=vfoptions; vfoptions_c.gridinterplayer=1; vfoptions_c.ngridinterp=5;
            simoptions_c=simoptions; simoptions_c.gridinterplayer=1; simoptions_c.ngridinterp=5;
            ngridinterp=5;
        end
        [~,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_zm,d_grid,a_grid,z_gridm,pi_zm,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions_c);
        StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_zm,pi_zm,simoptions_c,Params,[]);

        % The command, with unique-value TransitionProbs and then with quintile TransitionProbs
        simoptions_acp=simoptions_c;
        simoptions_acp.timehorizons=timehorizons;
        simoptions_acp.transprobs={'assets','shock'};
        ACP=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_zm,d_grid,a_grid,z_gridm,pi_zm,simoptions_acp);
        simoptions_q=simoptions_acp;
        simoptions_q.transprobs={'assets'};
        simoptions_q.transprobquantiles=nq;
        ACP_q=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_zm,d_grid,a_grid,z_gridm,pi_zm,simoptions_q);

        % Values of every function on the grid, (N_s,1), state s=(a,z,e) in column-major order
        ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy,FnsToEvaluate,Params,[],n_d,n_a,n_zm,d_grid,a_grid,z_gridm,simoptions_c);
        Vals=cell(length(FnNames),1);
        for ff=1:length(FnNames)
            Vals{ff}=gpuArray(reshape(gather(ValuesOnGrid.(FnNames{ff})),[N_s,1]));
        end
        dist=gpuArray(reshape(gather(StationaryDist),[N_s,1]));

        % One-period kernel T: T(s,s')=1{a'=aprime(s)} pi_z(z,z') pi_e(e') (two a' points with the GI weights)
        T=zeros(N_s,N_s);
        if n_semizm>0
            % T((a,semiz,z,e),(a',semiz',z',e'))=1{a'=aprime(a,semiz,z,e)} pi_semiz(semiz,semiz'|d2(a,semiz,z,e)) pi_z(z,z') pi_e(e')
            Policy_cpu=reshape(gather(Policy),[size(Policy,1),n_a,n_semizm,n_zr,n_er]); % rows (d2,aprime)
            for s1=1:n_semizm
                for zz=1:n_zr
                    for ee=1:n_er
                        for a1=1:n_a
                            row=a1+n_a*(s1-1)+n_a*n_semizm*(zz-1)+n_a*n_semizm*n_zr*(ee-1);
                            d2=Policy_cpu(1,a1,s1,zz,ee);
                            ap=Policy_cpu(2,a1,s1,zz,ee); % (the lower grid point, with the grid interpolation layer)
                            if ngridinterp==0
                                pup=0;
                            else
                                pup=(Policy_cpu(3,a1,s1,zz,ee)-1)/(ngridinterp+1); % rows (d2,aprime,L2,L2flag)
                                if Policy_cpu(4,a1,s1,zz,ee)==1
                                    pup=0; % all weight on the lower grid point
                                elseif Policy_cpu(4,a1,s1,zz,ee)==3
                                    pup=1; % all weight on the upper grid point
                                end
                            end
                            for s2=1:n_semizm
                                for zz2=1:n_zr
                                    for ee2=1:n_er
                                        prob=pi_semiz_bf(s1,s2,d2)*pi_z_cpu(zz,zz2);
                                        if n_em>0
                                            prob=prob*pi_e_cpu(ee2);
                                        end
                                        col=ap+n_a*(s2-1)+n_a*n_semizm*(zz2-1)+n_a*n_semizm*n_zr*(ee2-1);
                                        T(row,col)=T(row,col)+prob*(1-pup);
                                        if ngridinterp>0
                                            T(row,col+1)=T(row,col+1)+prob*pup; % upper point is ap+1, the next row of the same (semiz',z',e')
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        else
            Policy_cpu=reshape(gather(Policy),[size(Policy,1),n_a,n_zr,n_er]);
            for zz=1:n_zr
                for ee=1:n_er
                    rows=(1:n_a)'+n_a*(zz-1)+n_a*n_zr*(ee-1);
                    aplow=reshape(Policy_cpu(1,:,zz,ee),[n_a,1]);
                    if ngridinterp==0
                        pup=zeros(n_a,1);
                    else
                        L2index=reshape(Policy_cpu(2,:,zz,ee),[n_a,1]);
                        L2flag=reshape(Policy_cpu(3,:,zz,ee),[n_a,1]);
                        pup=(L2index-1)/(ngridinterp+1);
                        pup(L2flag==1)=0; % all weight on the lower grid point
                        pup(L2flag==3)=1; % all weight on the upper grid point
                    end
                    for zz2=1:n_zr
                        for ee2=1:n_er
                            prob=1;
                            if n_zm>0
                                prob=prob*pi_z_cpu(zz,zz2);
                            end
                            if n_em>0
                                prob=prob*pi_e_cpu(ee2);
                            end
                            cols=aplow+n_a*(zz2-1)+n_a*n_zr*(ee2-1);
                            T(rows+(cols-1)*N_s)=T(rows+(cols-1)*N_s)+prob*(1-pup);
                            if ngridinterp>0
                                T(rows+cols*N_s)=T(rows+cols*N_s)+prob*pup; % upper point is aplow+1, so one column further
                            end
                        end
                    end
                end
            end
        end
        T=gpuArray(T);
        % B1
        fprintf('AutoCorr brute force B1 (%s), kernel row sums minus one, this should be zero: %.3e \n',cfgstr,gather(max(abs(sum(T,2)-1))))
        fprintf('AutoCorr brute force B1 (%s), stationary dist pushed one period minus itself, this should be close to zero: %.3e \n',cfgstr,gather(max(abs(T'*dist-dist))))

        % B2. Moments
        Tk=eye(N_s,'gpuArray');
        for kk=1:Kmax
            Tk=Tk*T;
            hh=find(horizons==kk);
            if isempty(hh)
                continue
            end
            J=dist.*Tk; % joint distribution of (state now, state kk periods on), mass one
            for ff=1:length(FnNames)
                fn=FnNames{ff};
                x=Vals{ff};
                mu=sum(dist.*x);
                sd=sqrt(sum(dist.*(x-mu).^2));
                cov=sum(sum(J.*((x-mu)*(x-mu)'))); % the stationary dist is both marginals of J
                if kk==1
                    fprintf('AutoCorr brute force B2 (%s), %s Mean vs brute force, this should be zero: %.3e \n',cfgstr,fn,abs(gather(ACP.(fn).Mean)-gather(mu)))
                    fprintf('AutoCorr brute force B2 (%s), %s StdDeviation vs brute force, this should be zero: %.3e \n',cfgstr,fn,abs(gather(ACP.(fn).StdDeviation)-gather(sd)))
                    acov=gather(ACP.(fn).AutoCovariance); acorr=gather(ACP.(fn).AutoCorrelation);
                    hstr='';
                else
                    acov=gather(ACP.(fn).(['tperiods',num2str(kk)]).AutoCovariance); acorr=gather(ACP.(fn).(['tperiods',num2str(kk)]).AutoCorrelation);
                    hstr=sprintf(' tperiods%i',kk);
                end
                fprintf('AutoCorr brute force B2 (%s), %s%s AutoCovariance vs brute force, this should be zero: %.3e \n',cfgstr,fn,hstr,abs(acov-gather(cov)))
                fprintf('AutoCorr brute force B2 (%s), %s%s AutoCorrelation vs brute force, this should be zero: %.3e \n',cfgstr,fn,hstr,abs(acorr-gather(cov/sd^2)))
            end

            % B3. TransitionProbs between unique values (assets and shock)
            for fn={'assets','shock'}
                ff=find(strcmp(FnNames,fn{1}));
                [vv,~,idx]=unique(gather(Vals{ff}));
                nb=length(vv);
                S=sparse(1:N_s,idx,1,N_s,nb); % bin indicator
                S=gpuArray(full(S));
                massbin=S'*dist; % mass of each bin
                P_bf=gather((S'*J*S)./massbin); % row b: destination distribution of origin bin b (NaN where the bin has no mass)
                if kk==1
                    P=gather(ACP.(fn{1}).TransitionProbs);
                else
                    P=gather(ACP.(fn{1}).(['tperiods',num2str(kk)]).TransitionProbs);
                end
                if ~isequal(size(P),size(P_bf))
                    fprintf('AutoCorr brute force B3 (%s), %s%s TransitionProbs size [%s] differs from the brute force [%s], this should be zero: 1 \n',cfgstr,fn{1},hstr,num2str(size(P)),num2str(size(P_bf)))
                else
                    dd=abs(P-P_bf); dd=dd(~isnan(dd));
                    fprintf('AutoCorr brute force B3 (%s), %s%s TransitionProbs vs brute force, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,fn{1},hstr,max([dd(:);0]),sum(isnan(P(:))~=isnan(P_bf(:))))
                    rs=sum(P,2); rs=rs(gather(massbin)>0);
                    fprintf('AutoCorr brute force B3 (%s), %s%s TransitionProbs row sums minus one (bins with mass), this should be zero: %.3e \n',cfgstr,fn{1},hstr,max(abs(rs-1)))
                end
                % The shock's own transitions are known in closed form
                if strcmp(fn{1},'shock') && isequal(size(P),size(P_bf))
                    if mm==1 || mm==4 % z: the rows are those of pi_z^k (the unique values of z are the grid points, in the grid order)
                        P_an=pi_z_cpu^kk;
                        fprintf('AutoCorr brute force B3 (%s), shock%s TransitionProbs vs pi_z^k, this should be zero: %.3e \n',cfgstr,hstr,max(abs(P(:)-P_an(:))))
                    elseif mm==2 % e: every row is pi_e
                        P_an=repmat(pi_e_cpu',n_e,1);
                        fprintf('AutoCorr brute force B3 (%s), shock%s TransitionProbs vs pi_e in every row, this should be zero: %.3e \n',cfgstr,hstr,max(abs(P(:)-P_an(:))))
                    elseif nb==n_z*n_e % z*e: (z,e) to (z',e') is pi_z^k(z,z')*pi_e(e'), with the bins in the order of the sorted products
                        P_an=kron(pi_e_cpu',pi_z_cpu^kk); % rows (z,e) with z varying first, columns likewise
                        P_an=repmat(P_an,n_e,1); % the row only depends on z, so the kron's n_z rows repeat for every e
                        [~,ord]=sort(gather(Vals{ff}(1:n_a*n_semizr:N_s))); % the states (a=1,(semiz=1,)z,e) in column-major order (z first), sorted by the product: bin b is the pair ord(b)
                        P_an=P_an(ord,ord);
                        fprintf('AutoCorr brute force B3 (%s), shock%s TransitionProbs vs kron(pi_e,pi_z^k), this should be zero: %.3e \n',cfgstr,hstr,max(abs(P(:)-P_an(:))))
                    else
                        fprintf('AutoCorr brute force B3 (%s), shock%s: %i distinct products of z and e out of %i pairs, so the closed-form check is skipped \n',cfgstr,hstr,nb,n_z*n_e)
                    end
                end
            end

            % B4. TransitionProbs between quintiles (assets)
            ff=find(strcmp(FnNames,'assets'));
            x=gather(Vals{ff}); d_cpu=gather(dist);
            [xs,si]=sort(x);
            cm=cumsum(d_cpu(si));
            cutoffs=zeros(nq-1,1);
            for qq=1:nq-1
                cutoffs(qq)=xs(find(cm>qq/nq,1,'first'));
            end
            idx=ones(N_s,1);
            for qq=2:nq
                idx(x>cutoffs(qq-1))=qq;
            end
            S=gpuArray(full(sparse(1:N_s,idx,1,N_s,nq)));
            massbin=S'*dist;
            P_bf=gather((S'*J*S)./massbin);
            if kk==1
                P=gather(ACP_q.assets.TransitionProbs);
            else
                P=gather(ACP_q.assets.(['tperiods',num2str(kk)]).TransitionProbs);
            end
            dd=abs(P-P_bf); dd=dd(~isnan(dd));
            fprintf('AutoCorr brute force B4 (%s), assets%s quintile TransitionProbs vs brute force, this should be zero: %.3e, NaN-pattern mismatches %i \n',cfgstr,hstr,max([dd(:);0]),sum(isnan(P(:))~=isnan(P_bf(:))))
            fprintf('AutoCorr brute force B4 (%s), assets%s quintile TransitionProbs size [%i %i], expected [%i %i] \n',cfgstr,hstr,size(P,1),size(P,2),nq,nq)
        end
        clear T Tk J S
    end
end

output=struct();

end
