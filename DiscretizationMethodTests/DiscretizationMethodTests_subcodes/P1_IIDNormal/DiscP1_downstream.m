function output=DiscP1_downstream(calib)
% P1: feed this block's outputs into an actual toolkit solve
%
% Asserts only that the solve runs and that the objects come back the right shape. No economics,
% no statistics - the point is that SHAPE CONVENTIONS ARE ONLY WRONG RELATIVE TO A CONSUMER, and
% nothing else in this bank has a consumer. Three of the bugs found while scoping this bank were
% of exactly that kind: a life-cycle command returning no jequaloneDistz at all, two commands
% documenting the wrong grid shape, and MarkovChainMoments_FHorz demanding a transition slice it
% never read.
%
% An iid pi_e is fed in two ways here:
%   (a) disguised as a markov z with identical rows, which every InfHorz solver path accepts;
%   (b) NOT tested here: the toolkit's native iid e (vfoptions.n_e / e_grid / pi_e). The InfHorz
%       e path is out of scope for this bank - it belongs with whichever bank covers e - and the
%       CoreInfHorzTests e subcodes are themselves currently commented out.

fprintf('\n========== P1: downstream (do the outputs go into a solve) ========== \n')

output=struct();
mew=calib.mew; sigma=calib.sigma;
enum=9;

n_d=0; d_grid=[];
n_a=5; a_grid=linspace(0.1,5,n_a)'; % a few asset points here, unlike the n_a=1 model P0 uses
Params=struct(); Params.beta=0.9; Params.r=0.03;
DiscountFactorParamNames={'beta'};
ReturnFn=@(aprime,a,z,r) log(max((1+r)*a+exp(z)-aprime,10^(-10)));
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.zvalue=@(aprime,a,z) z;

cmdname={'IIDNormal_Tauchen','IIDNormal_TanakaToda','IID_Tauchen','IID_TanakaToda'};
for c_c=1:4
    if c_c==1
        [e_grid,pi_e]=discretizeIIDNormal_Tauchen(mew,sigma,enum,3,struct());
    elseif c_c==2
        [e_grid,pi_e]=discretizeIIDNormal_TanakaToda(mew,sigma,enum,struct());
    elseif c_c==3
        [e_grid,pi_e]=discretizeIID_Tauchen(mew,sigma,enum,3,struct());
    else
        [e_grid,pi_e]=discretizeIID_TanakaToda(mew,sigma,enum,struct());
    end

    % (a) as a markov z with identical rows
    z_grid=e_grid;
    pi_z=repmat(pi_e(:)',enum,1);
    n_z=enum;

    vfoptions=struct(); simoptions=struct();
    [V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params,[]);
    AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

    fprintf('%s: solve runs, size of V [T0], this should be zero: %i \n',cmdname{c_c},any(size(V)~=[n_a,n_z]))
    fprintf('%s: size of StationaryDist [T0], this should be zero: %i \n',cmdname{c_c},any(size(StationaryDist)~=[n_a,n_z]))
    fprintf('%s: StationaryDist sums to one [T1], this should be zero: %2.8e \n',cmdname{c_c},abs(sum(StationaryDist(:))-1))
    fprintf('%s: V is finite everywhere [T0], this should be zero: %i \n',cmdname{c_c},any(~isfinite(V(:))))

    % The z marginal of the stationary distribution must be pi_e itself: with identical rows the
    % chain is iid, so the distribution over z after one period IS pi_e regardless of assets.
    zmarginal=sum(StationaryDist,1);
    fprintf('%s: z marginal of the stationary dist equals pi_e [T1], this should be zero: %2.8e \n',cmdname{c_c},max(abs(zmarginal(:)-pi_e(:))))

    % And the mean of z from the toolkit's own statistics must match the chain's
    fprintf('%s: AllStats mean of z equals sum(pi_e.*e_grid) [T1], this should be zero: %2.8e \n',cmdname{c_c},abs(AllStats.zvalue.Mean-sum(pi_e(:).*e_grid(:))))

    output.cmd(c_c).name=cmdname{c_c};
    output.cmd(c_c).zmarginal=[zmarginal(:),pi_e(:)];
    output.cmd(c_c).meanz=[AllStats.zvalue.Mean,sum(pi_e(:).*e_grid(:))];
end

end
