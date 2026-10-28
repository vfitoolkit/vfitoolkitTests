function output=CoreInfHorzTPath_CrossTests_d1_semiz(T,PricePath,ParamPath,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline)
% Cross-tests for the InfHorz transition path with a semi-exogenous state, with d1 (labour supply) as well as d2.
% WRITTEN 2026-10-27, TEST-FIRST, against the InfHorz transition path semiz support being built in the
% toolkit. The transition path counterpart of CoreInfHorz_CrossTests2_d1_semiz in the stationary bank
% (CoreInfHorzTests): a semiz that is really just a markov -- its pi_semiz does not depend on d2 -- must
% give the same transition path as the same process declared as a markov z. Each runs the SAME path two
% ways that must agree:
%   (a) semiz as z: semiz with pi_semiz(:,:,d2)=pi_z for every d2 and semiz_grid=z_grid (no z) vs that
%       process as the markov z (no semiz)
%   (b) the same, with an additional z: (semiz,z) vs the two-dimensional markov z (zA,zB)
%   (c) the same as (a), with an iid e
%   (d) the same as (b), with an iid e
% and compares VPath, PolicyPath (no GI only), PolicyVals, AgentDistPath, the AggVars along the path, and
% one general eqm transition path (TransitionPath_InfHorz, Anderson, maxiter=3), each without and with
% grid interpolation. UNDER GI THE RAW POLICY INDEXES ARE NOT COMPARED, only the PolicyVals (two legal
% spellings of a point on a coarse node; see the doPart(1) subcodes).
%
% THE d2 ROW. The semiz side has a d2 the z side does not. Here d2 does nothing at all (pi_semiz does not
% depend on it, and searcheffortcost=0 takes it out of the ReturnFn), so it is an exact tie and the
% semiz side's d2 row is dropped before comparing: row 2 of the Policy and of the PolicyVals.
% The semiz side's Policy_final is handed to the z side with that row dropped too.
%
% pi_semiz IS GIVEN EXPLICITLY, as vfoptions.pi_semiz [n_semiz,n_semiz,n_d2] (the doPart(2) subcodes use
% the SemiExoStateFn route), so the two doParts between them cover both ways of declaring a semiz.
%
% BOTH SIDES OF EVERY COMPARISON ARE HANDED THE IDENTICAL V_final, Policy_final (less the d2 row) AND
% AgentDist_initial, solved once on the semiz side (the stationary semiz solve, without divide-and-
% conquer). In every disguise the two sides have the same state space in the same column-major order:
% semiz is the fastest exogenous state, which is also the first dimension of a two-dimensional z. So any
% difference is the transition path code's own. Agreement is to round-off, not bit-for-bit: the two code
% paths take the expectations in a different order.
%
% NOT IDENTICAL BY CONSTRUCTION: side 1 is the semiz code being built test-first, side 2 is the z code
% the transition path commands already support (with e in (c) and (d), the e code of doPart(1)/(8)).
% The additional z of (b) and (d) is a DIFFERENT process from the semiz one (3 points, rho=0.5, against
% the 5-point z), so a kron order swapped between the two sides changes the transition matrix and even
% the state count per dimension, and is caught.

Params.uempbenefit=0; % so the semiz ReturnFn is the z ReturnFn with semiz in place of z
Params.searcheffortcost=0; % d2 does nothing to the semiz transitions, and with no cost it does nothing at all
l_d1=1;
n_d1=n_d(1); d1_grid=d_grid(1:n_d1); % side 2 has d1 as its only decision variable
n_d2=n_d(end);

% The semiz process: the bank's markov z, declared as a semiz whose transitions ignore d2
n_semizX=n_z;
semiz_gridX=z_grid;
pi_semizX=repmat(pi_z,1,1,n_d2); % [n_semiz,n_semiz,n_d2], the same matrix for every d2

% The additional z of (b) and (d): a different process from the semiz one, so the kron order is detectable
n_zB=3;
[z_gridB,pi_zB]=discretizeAR1_FarmerToda(0,0.5,0.1,n_zB);
z_gridB=exp(z_gridB);

n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;

ReturnFn_semiz=@(d1,d2,aprime,a,semiz,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_noz_noe_semiz(d1,d2,aprime,a,semiz,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
ReturnFn_z=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
ReturnFn_semizz=@(d1,d2,aprime,a,semiz,z,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_z_noe_semiz(d1,d2,aprime,a,semiz,z,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
ReturnFn_zz=@(d,aprime,a,zA,zB,r,w,sigma,eta,varphi) ReturnFn_d_z_e_nosemiz(d,aprime,a,zA,zB,r,w,sigma,eta,varphi);
ReturnFn_semize=@(d1,d2,aprime,a,semiz,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_noz_e_semiz(d1,d2,aprime,a,semiz,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
ReturnFn_ze=@(d,aprime,a,z,e,r,w,sigma,eta,varphi) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,sigma,eta,varphi);
ReturnFn_semizze=@(d1,d2,aprime,a,semiz,z,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_z_e_semiz(d1,d2,aprime,a,semiz,z,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
ReturnFn_zze=@(d,aprime,a,zA,zB,e,r,w,sigma,eta,varphi) ReturnFn_d_z_e_nosemiz(d,aprime,a,zA*zB,e,r,w,sigma,eta,varphi);
% with uempbenefit=0 and searcheffortcost=0 each semiz ReturnFn is its z counterpart with semiz in place
% of the first z (all are c=(1+r)a+w*d*(product of the shocks)-aprime, plus the same leisure term)

FnsToEvaluate_semiz.assets=@(d1,d2,aprime,a,semiz) a;
FnsToEvaluate_semiz.earnings=@(d1,d2,aprime,a,semiz,w) w*d1*semiz;
FnsToEvaluate_semiz.nextassets=@(d1,d2,aprime,a,semiz) aprime;
FnsToEvaluate_semiz.shock=@(d1,d2,aprime,a,semiz) semiz;
FnsToEvaluate_z.assets=@(d,aprime,a,z) a;
FnsToEvaluate_z.earnings=@(d,aprime,a,z,w) w*d*z;
FnsToEvaluate_z.nextassets=@(d,aprime,a,z) aprime;
FnsToEvaluate_z.shock=@(d,aprime,a,z) z;
FnsToEvaluate_semizz.assets=@(d1,d2,aprime,a,semiz,z) a;
FnsToEvaluate_semizz.earnings=@(d1,d2,aprime,a,semiz,z,w) w*d1*semiz*z;
FnsToEvaluate_semizz.nextassets=@(d1,d2,aprime,a,semiz,z) aprime;
FnsToEvaluate_semizz.shock=@(d1,d2,aprime,a,semiz,z) semiz;
FnsToEvaluate_semizz.shockB=@(d1,d2,aprime,a,semiz,z) z;
FnsToEvaluate_zz.assets=@(d,aprime,a,zA,zB) a;
FnsToEvaluate_zz.earnings=@(d,aprime,a,zA,zB,w) w*d*zA*zB;
FnsToEvaluate_zz.nextassets=@(d,aprime,a,zA,zB) aprime;
FnsToEvaluate_zz.shock=@(d,aprime,a,zA,zB) zA;
FnsToEvaluate_zz.shockB=@(d,aprime,a,zA,zB) zB;
FnsToEvaluate_semize.assets=@(d1,d2,aprime,a,semiz,e) a;
FnsToEvaluate_semize.earnings=@(d1,d2,aprime,a,semiz,e,w) w*d1*semiz*e;
FnsToEvaluate_semize.nextassets=@(d1,d2,aprime,a,semiz,e) aprime;
FnsToEvaluate_semize.shock=@(d1,d2,aprime,a,semiz,e) semiz;
FnsToEvaluate_semize.eval=@(d1,d2,aprime,a,semiz,e) e;
FnsToEvaluate_ze.assets=@(d,aprime,a,z,e) a;
FnsToEvaluate_ze.earnings=@(d,aprime,a,z,e,w) w*d*z*e;
FnsToEvaluate_ze.nextassets=@(d,aprime,a,z,e) aprime;
FnsToEvaluate_ze.shock=@(d,aprime,a,z,e) z;
FnsToEvaluate_ze.eval=@(d,aprime,a,z,e) e;
FnsToEvaluate_semizze.assets=@(d1,d2,aprime,a,semiz,z,e) a;
FnsToEvaluate_semizze.earnings=@(d1,d2,aprime,a,semiz,z,e,w) w*d1*semiz*z*e;
FnsToEvaluate_semizze.nextassets=@(d1,d2,aprime,a,semiz,z,e) aprime;
FnsToEvaluate_semizze.shock=@(d1,d2,aprime,a,semiz,z,e) semiz;
FnsToEvaluate_semizze.shockB=@(d1,d2,aprime,a,semiz,z,e) z;
FnsToEvaluate_semizze.eval=@(d1,d2,aprime,a,semiz,z,e) e;
FnsToEvaluate_zze.assets=@(d,aprime,a,zA,zB,e) a;
FnsToEvaluate_zze.earnings=@(d,aprime,a,zA,zB,e,w) w*d*zA*zB*e;
FnsToEvaluate_zze.nextassets=@(d,aprime,a,zA,zB,e) aprime;
FnsToEvaluate_zze.shock=@(d,aprime,a,zA,zB,e) zA;
FnsToEvaluate_zze.shockB=@(d,aprime,a,zA,zB,e) zB;
FnsToEvaluate_zze.eval=@(d,aprime,a,zA,zB,e) e;

FnsToEvaluateGE_semiz.K=@(d1,d2,aprime,a,semiz) a;
FnsToEvaluateGE_semiz.N=@(d1,d2,aprime,a,semiz) d1*semiz;
FnsToEvaluateGE_z.K=@(d,aprime,a,z) a;
FnsToEvaluateGE_z.N=@(d,aprime,a,z) d*z;
FnsToEvaluateGE_semizz.K=@(d1,d2,aprime,a,semiz,z) a;
FnsToEvaluateGE_semizz.N=@(d1,d2,aprime,a,semiz,z) d1*semiz*z;
FnsToEvaluateGE_zz.K=@(d,aprime,a,zA,zB) a;
FnsToEvaluateGE_zz.N=@(d,aprime,a,zA,zB) d*zA*zB;
FnsToEvaluateGE_semize.K=@(d1,d2,aprime,a,semiz,e) a;
FnsToEvaluateGE_semize.N=@(d1,d2,aprime,a,semiz,e) d1*semiz*e;
FnsToEvaluateGE_ze.K=@(d,aprime,a,z,e) a;
FnsToEvaluateGE_ze.N=@(d,aprime,a,z,e) d*z*e;
FnsToEvaluateGE_semizze.K=@(d1,d2,aprime,a,semiz,z,e) a;
FnsToEvaluateGE_semizze.N=@(d1,d2,aprime,a,semiz,z,e) d1*semiz*z*e;
FnsToEvaluateGE_zze.K=@(d,aprime,a,zA,zB,e) a;
FnsToEvaluateGE_zze.N=@(d,aprime,a,zA,zB,e) d*zA*zB*e;

% General eqm: the same Cobb-Douglas firm and Anderson setup as the doPart(3) cross-tests. The path starts
% from a flat guess near the general eqm (see below), which is NOT an equilibrium path, and is cut off at
% maxiter=3: the point is to have the two sides run the same three iterations, so that every update -- the
% agent dist, the aggregates, the GE eqns and the Anderson step -- has to agree, not just the first.
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;
GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));
% THE INITIAL GUESS IS NEAR THE GENERAL EQM, not the partial eqm prices of the rest of the bank (w=1): from w=1 the
% wage residual is ~0.56, Anderson extrapolates on it, and on 2026-10-28 it took test E(z) to NaN prices at
% iteration 4. r=0.052, w=0.44 is about where the z models' stationary general eqm sits (doPart(1)). The test
% only needs both sides to take the same steps, and the prices to move (checked below), not a converged path,
% so maxiter is also kept short.
PricePathGE.r=0.052*ones(1,T);
PricePathGE.w=0.44*ones(1,T);
ParamPathGE.sigma=Params.sigma*ones(1,T); % w is a general eqm price here, so it cannot also be on the ParamPath
transpathoptionsGE=transpathoptionsbaseline;
transpathoptionsGE.GEnewprice=2;
transpathoptionsGE.GEnewprice2.howtoupdate={'CapitalMarket','r',0,0.1; ...
                                            'LabourMarket','w',0,0.1};
transpathoptionsGE.maxiter=3;
transpathoptionsGE.verbose=0;

for cc=1:4
    % side 1 always has the semiz (via vfoptions/simoptions, with pi_semiz given explicitly); side 2 is the
    % same model with that process as (the first dimension of) a markov z
    vfo1=struct(); simo1=struct(); vfo2=struct(); simo2=struct();
    vfo1.n_semiz=n_semizX; vfo1.semiz_grid=semiz_gridX; vfo1.pi_semiz=pi_semizX; vfo1.l_dsemiz=1;
    simo1.n_semiz=n_semizX; simo1.semiz_grid=semiz_gridX; simo1.pi_semiz=pi_semizX; simo1.l_dsemiz=1;
    simo1.d_grid=d_grid;
    switch cc
        case 1
            testname='Cross test S(a) (semiz as z)';
            n_z1=0; z_grid1=[]; pi_z1=[];
            ReturnFn1=ReturnFn_semiz; FnsToEvaluate1=FnsToEvaluate_semiz; FnsToEvaluateGE1=FnsToEvaluateGE_semiz;
            n_z2=n_semizX; z_grid2=semiz_gridX; pi_z2=pi_z;
            ReturnFn2=ReturnFn_z; FnsToEvaluate2=FnsToEvaluate_z; FnsToEvaluateGE2=FnsToEvaluateGE_z;
        case 2
            testname='Cross test S(b) (semiz as z, with z)';
            n_z1=n_zB; z_grid1=z_gridB; pi_z1=pi_zB;
            ReturnFn1=ReturnFn_semizz; FnsToEvaluate1=FnsToEvaluate_semizz; FnsToEvaluateGE1=FnsToEvaluateGE_semizz;
            % (semiz,z) as one two-dimensional markov z, n_z=[n_semiz,n_zB], the semiz process first (fastest),
            % as semiz is in the state order (semiz,z). The joint index is s+n_semiz*(zB-1), the same order as
            % the (semiz,z) states of side 1. Its transition matrix has entry ((s,zB),(s',zB')) =
            % pi_semiz(s,s')*pi_zB(zB,zB'), which is kron(pi_zB,pi_z): kron(A,B)(k+n_B*(i-1), l+n_B*(j-1)) =
            % A(i,j)*B(k,l), so with A=pi_zB (index i,j over zB) and B=pi_z (index k,l over semiz, n_B=n_semiz)
            % row (s=k,zB=i), column (s'=l,zB'=j) is pi_zB(i,j)*pi_z(k,l). The grid is stacked, semiz on top.
            % kron(pi_z,pi_zB) would be the WRONG order; pi_zB differs from pi_z in both values and size, so
            % that mistake changes the problem and is caught.
            n_z2=[n_semizX,n_zB]; z_grid2=[semiz_gridX;z_gridB]; pi_z2=kron(pi_zB,pi_z);
            ReturnFn2=ReturnFn_zz; FnsToEvaluate2=FnsToEvaluate_zz; FnsToEvaluateGE2=FnsToEvaluateGE_zz;
        case 3
            testname='Cross test S(c) (semiz as z, with e)';
            n_z1=0; z_grid1=[]; pi_z1=[];
            vfo1.n_e=n_e; vfo1.e_grid=e_grid; vfo1.pi_e=pi_e;
            simo1.n_e=n_e; simo1.e_grid=e_grid; simo1.pi_e=pi_e;
            ReturnFn1=ReturnFn_semize; FnsToEvaluate1=FnsToEvaluate_semize; FnsToEvaluateGE1=FnsToEvaluateGE_semize;
            n_z2=n_semizX; z_grid2=semiz_gridX; pi_z2=pi_z;
            vfo2.n_e=n_e; vfo2.e_grid=e_grid; vfo2.pi_e=pi_e;
            simo2.n_e=n_e; simo2.e_grid=e_grid; simo2.pi_e=pi_e;
            ReturnFn2=ReturnFn_ze; FnsToEvaluate2=FnsToEvaluate_ze; FnsToEvaluateGE2=FnsToEvaluateGE_ze;
        case 4
            testname='Cross test S(d) (semiz as z, with z and e)';
            n_z1=n_zB; z_grid1=z_gridB; pi_z1=pi_zB;
            vfo1.n_e=n_e; vfo1.e_grid=e_grid; vfo1.pi_e=pi_e;
            simo1.n_e=n_e; simo1.e_grid=e_grid; simo1.pi_e=pi_e;
            ReturnFn1=ReturnFn_semizze; FnsToEvaluate1=FnsToEvaluate_semizze; FnsToEvaluateGE1=FnsToEvaluateGE_semizze;
            % the same two-dimensional z as S(b) (kron order derived there); e is unchanged, after both
            n_z2=[n_semizX,n_zB]; z_grid2=[semiz_gridX;z_gridB]; pi_z2=kron(pi_zB,pi_z);
            vfo2.n_e=n_e; vfo2.e_grid=e_grid; vfo2.pi_e=pi_e;
            simo2.n_e=n_e; simo2.e_grid=e_grid; simo2.pi_e=pi_e;
            ReturnFn2=ReturnFn_zze; FnsToEvaluate2=FnsToEvaluate_zze; FnsToEvaluateGE2=FnsToEvaluateGE_zze;
    end

    for gi=0:1
        vfoptions1=vfo1; simoptions1=simo1; vfoptions2=vfo2; simoptions2=simo2;
        if gi==0
            tiername='no GI';
        else
            tiername='with GI';
            vfoptions1.gridinterplayer=1; vfoptions1.ngridinterp=5;
            simoptions1.gridinterplayer=1; simoptions1.ngridinterp=5;
            vfoptions2.gridinterplayer=1; vfoptions2.ngridinterp=5;
            simoptions2.gridinterplayer=1; simoptions2.ngridinterp=5;
        end

        % V_final, Policy_final and AgentDist_initial, solved once on side 1 (stationary semiz, no DC) and
        % handed to both sides; side 2 gets Policy_final without the d2 row
        [V_final,Policy_final]=ValueFnIter_InfHorz(n_d,n_a,n_z1,d_grid,a_grid,z_grid1,pi_z1,ReturnFn1,Params,DiscountFactorParamNames,[],vfoptions1);
        AgentDist_initial=StationaryDist_InfHorz(Policy_final,n_d,n_a,n_z1,pi_z1,simoptions1,Params,[]);
        szP=size(Policy_final);
        keeprows=setdiff(1:szP(1),l_d1+1); % every row but d2
        Policy_final2=reshape(Policy_final(keeprows,:),[numel(keeprows),szP(2:end)]);

        [VPath1,PolicyPath1]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z1, d_grid, a_grid,z_grid1, pi_z1, DiscountFactorParamNames, ReturnFn1, transpathoptionsbaseline, vfoptions1);
        [VPath2,PolicyPath2]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final2, Params, n_d1, n_a, n_z2, d1_grid, a_grid,z_grid2, pi_z2, DiscountFactorParamNames, ReturnFn2, transpathoptionsbaseline, vfoptions2);
        fprintf('%s, %s, VPath, this should be zero: %.3e \n',testname,tiername,max(abs(gather(VPath1(:))-gather(VPath2(:)))))
        if gi==0
            PolicyPath1short=PolicyPath1(keeprows,:); % drop the d2 row (rows are the first dim, the rest in the same order)
            fprintf('%s, %s, PolicyPath (less the d2 row), this should be zero: %.3e \n',testname,tiername,max(abs(gather(PolicyPath1short(:))-gather(PolicyPath2(:)))))
            clear PolicyPath1short
        end
        PolicyValsPath1=PolicyInd2Val_InfHorz_TPath(PolicyPath1,n_d,n_a,n_z1,T,d_grid,a_grid,vfoptions1);
        PolicyValsPath2=PolicyInd2Val_InfHorz_TPath(PolicyPath2,n_d1,n_a,n_z2,T,d1_grid,a_grid,vfoptions2);
        szPV=size(PolicyValsPath1);
        PolicyValsPath1short=PolicyValsPath1(setdiff(1:szPV(1),l_d1+1),:); % drop the d2 value row
        fprintf('%s, %s, PolicyVals along the path (less the d2 row), this should be zero: %.3e \n',testname,tiername,max(abs(gather(PolicyValsPath1short(:))-gather(PolicyValsPath2(:)))))

        AgentDistPath1=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPath, PolicyPath1, n_d, n_a, n_z1, pi_z1, T, Params, transpathoptionsbaseline, simoptions1);
        AgentDistPath2=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPath, PolicyPath2, n_d1, n_a, n_z2, pi_z2, T, Params, transpathoptionsbaseline, simoptions2);
        fprintf('%s, %s, AgentDistPath, this should be zero: %.3e \n',testname,tiername,max(abs(gather(AgentDistPath1(:))-gather(AgentDistPath2(:)))))

        AggVarsPath1=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate1, AgentDistPath1, PolicyPath1, PricePath, ParamPath, Params, T, n_d, n_a, n_z1, d_grid, a_grid, z_grid1, transpathoptionsbaseline, simoptions1);
        AggVarsPath2=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate2, AgentDistPath2, PolicyPath2, PricePath, ParamPath, Params, T, n_d1, n_a, n_z2, d1_grid, a_grid, z_grid2, transpathoptionsbaseline, simoptions2);
        for fnname=fieldnames(FnsToEvaluate2)'
            fprintf('%s, %s, AggVars %s Mean along the path, this should be zero: %.3e \n',testname,tiername,fnname{1},max(abs(gather(AggVarsPath1.(fnname{1}).Mean(:))-gather(AggVarsPath2.(fnname{1}).Mean(:)))))
        end

        [PricePathOut1,GEcondnPath1]=TransitionPath_InfHorz(PricePathGE, ParamPathGE, T, V_final, AgentDist_initial, n_d, n_a, n_z1, d_grid,a_grid,z_grid1, pi_z1, ReturnFn1, FnsToEvaluateGE1, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, transpathoptionsGE, simoptions1, vfoptions1);
        [PricePathOut2,GEcondnPath2]=TransitionPath_InfHorz(PricePathGE, ParamPathGE, T, V_final, AgentDist_initial, n_d1, n_a, n_z2, d1_grid,a_grid,z_grid2, pi_z2, ReturnFn2, FnsToEvaluateGE2, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, transpathoptionsGE, simoptions2, vfoptions2);
        fprintf('%s, %s, general eqm path (Anderson, maxiter=3), this should be zero, PricePath r: %.3e, w: %.3e \n',testname,tiername,max(abs(gather(PricePathOut1.r(:))-gather(PricePathOut2.r(:)))),max(abs(gather(PricePathOut1.w(:))-gather(PricePathOut2.w(:)))))
        fprintf('%s, %s, general eqm path (Anderson, maxiter=3), this should be zero, GEcondnPath CapitalMarket: %.3e, LabourMarket: %.3e \n',testname,tiername,max(abs(gather(GEcondnPath1.CapitalMarket(:))-gather(GEcondnPath2.CapitalMarket(:)))),max(abs(gather(GEcondnPath1.LabourMarket(:))-gather(GEcondnPath2.LabourMarket(:)))))
        % Non-vacuity: the five iterations must actually have moved the prices, else agreeing on them says little
        fprintf('%s, %s, general eqm path moved the prices away from the initial guess, this should NOT be zero: %.3e \n',testname,tiername,max(abs(gather(PricePathOut1.r(:))-PricePathGE.r(:))))

        clear V_final Policy_final Policy_final2 AgentDist_initial szP keeprows szPV VPath1 VPath2 PolicyPath1 PolicyPath2 PolicyValsPath1 PolicyValsPath1short PolicyValsPath2 AgentDistPath1 AgentDistPath2 AggVarsPath1 AggVarsPath2 PricePathOut1 PricePathOut2 GEcondnPath1 GEcondnPath2
    end
end

%%
output=struct();

end
