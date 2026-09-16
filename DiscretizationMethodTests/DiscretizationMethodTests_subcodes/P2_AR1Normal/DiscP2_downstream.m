function output=DiscP2_downstream(calib)
% P2: feed this block's outputs into an actual toolkit solve
%
% Asserts only that the solve runs and that the objects come back the right shape, plus one
% substantive check that a shape-only test cannot give: the z marginal of the stationary
% distribution must equal the chain's own stationary distribution. If a grid or transition matrix
% were being consumed in the wrong orientation, the solve would still run and still produce a
% distribution - it would just be the wrong one.

fprintf('\n========== P2: downstream (do the outputs go into a solve) ========== \n')

output=struct();
znum=9;
n_d=0; d_grid=[];
n_a=5; a_grid=linspace(0.1,5,n_a)';
Params=struct(); Params.beta=0.9; Params.r=0.03;
DiscountFactorParamNames={'beta'};
ReturnFn=@(aprime,a,z,r) log(max((1+r)*a+exp(z)-aprime,10^(-10)));
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.zvalue=@(aprime,a,z) z;

cmdname={'AR1_Tauchen','AR1_Rouwenhorst','AR1_TauchenHussey','AR1_FarmerToda'};
for cal_c=1:3
    cname=calib.names{cal_c};
    mew=calib.(cname).mew; rho=calib.(cname).rho; sigma=calib.(cname).sigma;
    fprintf('\n--- calibration %s --- \n',cname)
    for c_c=1:4
        if c_c==1
            [z_grid,pi_z]=discretizeAR1_Tauchen(mew,rho,sigma,znum,3,struct());
        elseif c_c==2
            [z_grid,pi_z]=discretizeAR1_Rouwenhorst(mew,rho,sigma,znum,struct());
        elseif c_c==3
            [z_grid,pi_z]=discretizeAR1_TauchenHussey(mew,rho,sigma,znum,struct());
        else
            [z_grid,pi_z,~]=discretizeAR1_FarmerToda(mew,rho,sigma,znum,struct('verbose',0));
        end
        n_z=znum;

        vfoptions=struct(); simoptions=struct();
        [V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
        StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params,[]);
        AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

        fprintf('%s/%s: solve runs, size of V [T0], this should be zero: %i \n',cname,cmdname{c_c},any(size(V)~=[n_a,n_z]))
        fprintf('%s/%s: StationaryDist sums to one [T1], this should be zero: %2.8e \n',cname,cmdname{c_c},abs(sum(StationaryDist(:))-1))
        fprintf('%s/%s: V is finite everywhere [T0], this should be zero: %i \n',cname,cmdname{c_c},any(~isfinite(V(:))))

        % The substantive one: the z marginal of the agent distribution must be the chain's own
        % stationary distribution. The asset dimension cannot affect it, because z is exogenous.
        [~,~,~,mcstatdist]=MarkovChainMoments(z_grid,pi_z);
        zmarg=sum(StationaryDist,1);
        fprintf('%s/%s: z marginal of the agent dist equals the chain stationary dist [T2], this should be small: %2.8e \n',cname,cmdname{c_c},max(abs(zmarg(:)-mcstatdist(:))))
        % T2 rather than T1: StationaryDist_InfHorz iterates to its own convergence tolerance, so
        % the agreement is to that tolerance rather than to machine precision.
        fprintf('%s/%s: AllStats mean of z equals the chain mean [T2], this should be small: %2.8e \n',cname,cmdname{c_c},abs(AllStats.zvalue.Mean-sum(mcstatdist(:).*z_grid(:))))

        output.(cname).cmd(c_c).name=cmdname{c_c};
        output.(cname).cmd(c_c).zmarginalerr=max(abs(zmarg(:)-mcstatdist(:)));
    end
end

end
