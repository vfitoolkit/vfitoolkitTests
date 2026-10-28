function output=CoreInfHorzTPath_CrossTests_nod_nosemiz_e(T,PricePath,ParamPath,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline)
% Cross-tests for the InfHorz transition path with an iid e, without d: D, E and G.
% The transition path counterpart of CoreInfHorz_CrossTests_nod_nosemiz_e in the stationary bank
% (CoreInfHorzTests), with the same letters for the same disguises. Each runs the SAME transition
% path two ways that must agree:
%   (D) e-as-z: an iid e (no z) vs the very same process declared as a markov z whose every row is pi_e'
%   (E) trivial e: n_e=1, e_grid=1, pi_e=1 vs the same model with no e. Done twice, on the noz model
%       and on the z model, so that the trivial e is checked both alone and alongside a markov z.
%   (G) (z,e) vs the same two shocks merged into one two-dimensional markov z
% and compares VPath, PolicyPath, AgentDistPath, the AggVars along the path, and one general eqm
% transition path (TransitionPath_InfHorz, Anderson). Each without and with grid interpolation
% (gridinterplayer=1, ngridinterp=5). UNDER GI THE RAW POLICY INDEXES ARE NOT COMPARED, only the
% PolicyVals: a point on a coarse node has two legal spellings under the interpolation layer, so two
% correct implementations need not return the same index (see the doPart(1) subcodes).
% There is no (F), the e-permutation test of the stationary bank: it is a statement about relabelling,
% which that bank covers, and nothing in it is specific to the path.
%
% BOTH SIDES OF EVERY COMPARISON ARE HANDED THE IDENTICAL V_final, Policy_final AND AgentDist_initial.
% They are solved once, on side 1, and passed to side 2 as the same array. In every disguise the two
% sides have the same state space in the same column-major order -- a, then z, then e, which is also
% the order of a two-dimensional z with n_z=[n_z,n_e] -- so the arrays need no reshape. That way any
% difference is the transition path code's own, and not the stationary solver's (the stationary bank
% tests that).
%
% Agreement is to round-off, not bit-for-bit: the two code paths take the expectation over the shocks
% in a different order. A raw Policy index could in principle flip on an exact tie at round-off; the
% VPath line beside it says whether that is what happened.
%
% Side 1 is the e model in every test, i.e. the code being built test-first; side 2 is a model the
% transition path commands already support (z without e), except in E(noz), where side 2 is the noz
% noe model, which is itself new.

n_d=0; d_grid=[];
n_z_withz=n_z; z_grid_withz=z_grid; pi_z_withz=pi_z; % for E(z) and G

n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;

% All four are the same formula, c=(1+r)a+w*(product of the shocks)-aprime, so at the same shock grid
% and the same shock process they give identical F
ReturnFn_none=@(aprime,a,r,w,sigma) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,sigma);
ReturnFn_e=@(aprime,a,e,r,w,sigma) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,sigma);
ReturnFn_z=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
ReturnFn_ze=@(aprime,a,z,e,r,w,sigma) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,sigma);

FnsToEvaluate_none.assets=@(aprime,a) a;
FnsToEvaluate_none.earnings=@(aprime,a,w) w;
FnsToEvaluate_none.nextassets=@(aprime,a) aprime;
FnsToEvaluate_e.assets=@(aprime,a,e) a;
FnsToEvaluate_e.earnings=@(aprime,a,e,w) w*e;
FnsToEvaluate_e.nextassets=@(aprime,a,e) aprime;
FnsToEvaluate_e.shock=@(aprime,a,e) e;
FnsToEvaluate_z.assets=@(aprime,a,z) a;
FnsToEvaluate_z.earnings=@(aprime,a,z,w) w*z;
FnsToEvaluate_z.nextassets=@(aprime,a,z) aprime;
FnsToEvaluate_z.shock=@(aprime,a,z) z;
FnsToEvaluate_ze.assets=@(aprime,a,z,e) a;
FnsToEvaluate_ze.earnings=@(aprime,a,z,e,w) w*z*e;
FnsToEvaluate_ze.nextassets=@(aprime,a,z,e) aprime;
FnsToEvaluate_ze.zval=@(aprime,a,z,e) z;
FnsToEvaluate_ze.eval=@(aprime,a,z,e) e; % in G, side 2 hands the second z in this position
% Test E without z: the trivial-e side takes e as an input (the toolkit passes it, whatever its value), so it needs
% its own FnsToEvaluate with an e. Same three fields as FnsToEvaluate_none, and with e=1 the same values.
FnsToEvaluate_etriv.assets=@(aprime,a,e) a;
FnsToEvaluate_etriv.earnings=@(aprime,a,e,w) w*e;
FnsToEvaluate_etriv.nextassets=@(aprime,a,e) aprime;
% E(z): the z and trivial e model, with the fields of FnsToEvaluate_z so the two sides can be compared field by field
FnsToEvaluate_zetriv.assets=@(aprime,a,z,e) a;
FnsToEvaluate_zetriv.earnings=@(aprime,a,z,e,w) w*z*e;
FnsToEvaluate_zetriv.nextassets=@(aprime,a,z,e) aprime;
FnsToEvaluate_zetriv.shock=@(aprime,a,z,e) z;

% General eqm: the same Cobb-Douglas firm as the doPart(1) subcodes, both prices from its first order
% conditions, and the same Anderson setup (GEnewprice=2, residual-form eqns, howtoupdate rows in the
% order of the PricePath fields). The path starts from the doPart(1) PricePath.r and ParamPath.w, which
% is NOT an equilibrium path, and is cut off at maxiter: the point is not to solve the general eqm but
% to have the two sides run the same five iterations, so that every update -- the agent dist, the
% aggregates, the GE eqns and the Anderson step -- has to agree, not just the first.
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;
GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));
FnsToEvaluateGE_e.K=@(aprime,a,e) a;
FnsToEvaluateGE_e.N=@(aprime,a,e) e;
FnsToEvaluateGE_z.K=@(aprime,a,z) a;
FnsToEvaluateGE_z.N=@(aprime,a,z) z;
FnsToEvaluateGE_ze.K=@(aprime,a,z,e) a;
FnsToEvaluateGE_ze.N=@(aprime,a,z,e) z*e;
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
    % side 1 always has the iid e; side 2 is the same model in disguise
    vfo1=struct(); simo1=struct(); vfo2=struct(); simo2=struct();
    switch cc
        case 1
            testname='Cross test D (iid e == same iid as markov z)';
            n_z1=0; z_grid1=[]; pi_z1=[];
            vfo1.n_e=n_e; vfo1.e_grid=e_grid; vfo1.pi_e=pi_e;
            simo1.n_e=n_e; simo1.e_grid=e_grid; simo1.pi_e=pi_e;
            ReturnFn1=ReturnFn_e; FnsToEvaluate1=FnsToEvaluate_e; FnsToEvaluateGE1=FnsToEvaluateGE_e;
            % The same iid shock declared as a markov z: every row of pi_z is pi_e', so z is iid with distribution pi_e
            n_z2=n_e; z_grid2=e_grid; pi_z2=repmat(pi_e',n_e,1);
            ReturnFn2=ReturnFn_z; FnsToEvaluate2=FnsToEvaluate_z; FnsToEvaluateGE2=FnsToEvaluateGE_z;
            doGE=1;
        case 2
            testname='Cross test E, no z (trivial e == noe)';
            n_z1=0; z_grid1=[]; pi_z1=[];
            vfo1.n_e=1; vfo1.e_grid=1; vfo1.pi_e=1;
            simo1.n_e=1; simo1.e_grid=1; simo1.pi_e=1;
            ReturnFn1=ReturnFn_e; FnsToEvaluate1=FnsToEvaluate_etriv; FnsToEvaluateGE1=struct();
            n_z2=0; z_grid2=[]; pi_z2=[];
            ReturnFn2=ReturnFn_none; FnsToEvaluate2=FnsToEvaluate_none; FnsToEvaluateGE2=struct();
            % NO GENERAL EQM IN THIS ONE. With no shock at all and beta*(1+r)<1 everyone runs their assets
            % down to a=0, so K=0 and the firm's MPK is infinite: the GE eqns are not defined on this path.
            % The trivial e meets the general eqm in E(z) below instead.
            doGE=0;
        case 3
            testname='Cross test E, with z (z and trivial e == z noe)';
            n_z1=n_z_withz; z_grid1=z_grid_withz; pi_z1=pi_z_withz;
            vfo1.n_e=1; vfo1.e_grid=1; vfo1.pi_e=1;
            simo1.n_e=1; simo1.e_grid=1; simo1.pi_e=1;
            ReturnFn1=ReturnFn_ze; FnsToEvaluate1=FnsToEvaluate_zetriv; FnsToEvaluateGE1=FnsToEvaluateGE_ze; % with e=1, N=z*e is N=z
            n_z2=n_z_withz; z_grid2=z_grid_withz; pi_z2=pi_z_withz;
            ReturnFn2=ReturnFn_z; FnsToEvaluate2=FnsToEvaluate_z; FnsToEvaluateGE2=FnsToEvaluateGE_z;
            doGE=1;
        case 4
            testname='Cross test G ((z,e) == two-dim markov z)';
            n_z1=n_z_withz; z_grid1=z_grid_withz; pi_z1=pi_z_withz;
            vfo1.n_e=n_e; vfo1.e_grid=e_grid; vfo1.pi_e=pi_e;
            simo1.n_e=n_e; simo1.e_grid=e_grid; simo1.pi_e=pi_e;
            ReturnFn1=ReturnFn_ze; FnsToEvaluate1=FnsToEvaluate_ze; FnsToEvaluateGE1=FnsToEvaluateGE_ze;
            % The same two shocks as one two-dimensional markov z, n_z=[n_z,n_e], z varying fastest: the
            % joint index is z+n_z*(e-1), the same order as the (z,e) states of side 1. Its transition matrix
            % has entry ((z,e),(z',e')) = pi_z(z,z')*pi_e(e'), which is kron(A,pi_z) with A=ones(n_e,1)*pi_e',
            % since kron(A,B)(k+n_z*(i-1), l+n_z*(j-1)) = A(i,j)*B(k,l) = pi_e(j)*pi_z(k,l) for row (z=k,e=i)
            % and column (z'=l,e'=j). The rows of A are all pi_e', so e' does not depend on e: iid, as it must.
            % The grid is stacked, z_grid on top of e_grid. The ReturnFn and FnsToEvaluate are unchanged,
            % the e input is just the second z.
            n_z2=[n_z_withz,n_e]; z_grid2=[z_grid_withz;e_grid]; pi_z2=kron(ones(n_e,1)*pi_e',pi_z_withz);
            ReturnFn2=ReturnFn_ze; FnsToEvaluate2=FnsToEvaluate_ze; FnsToEvaluateGE2=FnsToEvaluateGE_ze;
            doGE=1;
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

        % V_final, Policy_final and AgentDist_initial, solved once on side 1 and handed to both sides
        [V_final,Policy_final]=ValueFnIter_InfHorz(n_d,n_a,n_z1,d_grid,a_grid,z_grid1,pi_z1,ReturnFn1,Params,DiscountFactorParamNames,[],vfoptions1);
        AgentDist_initial=StationaryDist_InfHorz(Policy_final,n_d,n_a,n_z1,pi_z1,simoptions1,Params,[]);

        [VPath1,PolicyPath1]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z1, d_grid, a_grid,z_grid1, pi_z1, DiscountFactorParamNames, ReturnFn1, transpathoptionsbaseline, vfoptions1);
        [VPath2,PolicyPath2]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z2, d_grid, a_grid,z_grid2, pi_z2, DiscountFactorParamNames, ReturnFn2, transpathoptionsbaseline, vfoptions2);
        fprintf('%s, %s, VPath, this should be zero: %.3e \n',testname,tiername,max(abs(gather(VPath1(:))-gather(VPath2(:)))))
        if gi==0
            fprintf('%s, %s, PolicyPath, this should be zero: %.3e \n',testname,tiername,max(abs(gather(PolicyPath1(:))-gather(PolicyPath2(:)))))
        end
        PolicyValsPath1=PolicyInd2Val_InfHorz_TPath(PolicyPath1,n_d,n_a,n_z1,T,d_grid,a_grid,vfoptions1);
        PolicyValsPath2=PolicyInd2Val_InfHorz_TPath(PolicyPath2,n_d,n_a,n_z2,T,d_grid,a_grid,vfoptions2);
        fprintf('%s, %s, PolicyVals along the path, this should be zero: %.3e \n',testname,tiername,max(abs(gather(PolicyValsPath1(:))-gather(PolicyValsPath2(:)))))

        AgentDistPath1=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPath, PolicyPath1, n_d, n_a, n_z1, pi_z1, T, Params, transpathoptionsbaseline, simoptions1);
        AgentDistPath2=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPath, PolicyPath2, n_d, n_a, n_z2, pi_z2, T, Params, transpathoptionsbaseline, simoptions2);
        fprintf('%s, %s, AgentDistPath, this should be zero: %.3e \n',testname,tiername,max(abs(gather(AgentDistPath1(:))-gather(AgentDistPath2(:)))))

        AggVarsPath1=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate1, AgentDistPath1, PolicyPath1, PricePath, ParamPath, Params, T, n_d, n_a, n_z1, d_grid, a_grid, z_grid1, transpathoptionsbaseline, simoptions1);
        AggVarsPath2=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate2, AgentDistPath2, PolicyPath2, PricePath, ParamPath, Params, T, n_d, n_a, n_z2, d_grid, a_grid, z_grid2, transpathoptionsbaseline, simoptions2);
        for fnname=fieldnames(FnsToEvaluate2)'
            fprintf('%s, %s, AggVars %s Mean along the path, this should be zero: %.3e \n',testname,tiername,fnname{1},max(abs(gather(AggVarsPath1.(fnname{1}).Mean(:))-gather(AggVarsPath2.(fnname{1}).Mean(:)))))
        end

        if doGE==1
            [PricePathOut1,GEcondnPath1]=TransitionPath_InfHorz(PricePathGE, ParamPathGE, T, V_final, AgentDist_initial, n_d, n_a, n_z1, d_grid,a_grid,z_grid1, pi_z1, ReturnFn1, FnsToEvaluateGE1, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, transpathoptionsGE, simoptions1, vfoptions1);
            [PricePathOut2,GEcondnPath2]=TransitionPath_InfHorz(PricePathGE, ParamPathGE, T, V_final, AgentDist_initial, n_d, n_a, n_z2, d_grid,a_grid,z_grid2, pi_z2, ReturnFn2, FnsToEvaluateGE2, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, transpathoptionsGE, simoptions2, vfoptions2);
            fprintf('%s, %s, general eqm path (Anderson, maxiter=3), this should be zero, PricePath r: %.3e, w: %.3e \n',testname,tiername,max(abs(gather(PricePathOut1.r(:))-gather(PricePathOut2.r(:)))),max(abs(gather(PricePathOut1.w(:))-gather(PricePathOut2.w(:)))))
            fprintf('%s, %s, general eqm path (Anderson, maxiter=3), this should be zero, GEcondnPath CapitalMarket: %.3e, LabourMarket: %.3e \n',testname,tiername,max(abs(gather(GEcondnPath1.CapitalMarket(:))-gather(GEcondnPath2.CapitalMarket(:)))),max(abs(gather(GEcondnPath1.LabourMarket(:))-gather(GEcondnPath2.LabourMarket(:)))))
            % Non-vacuity: the five iterations must actually have moved the prices, else agreeing on them says little
            fprintf('%s, %s, general eqm path moved the prices away from the initial guess, this should NOT be zero: %.3e \n',testname,tiername,max(abs(gather(PricePathOut1.r(:))-PricePathGE.r(:))))
            clear PricePathOut1 PricePathOut2 GEcondnPath1 GEcondnPath2
        end

        clear V_final Policy_final AgentDist_initial VPath1 VPath2 PolicyPath1 PolicyPath2 PolicyValsPath1 PolicyValsPath2 AgentDistPath1 AgentDistPath2 AggVarsPath1 AggVarsPath2
    end
end

%%
output=struct();

end
