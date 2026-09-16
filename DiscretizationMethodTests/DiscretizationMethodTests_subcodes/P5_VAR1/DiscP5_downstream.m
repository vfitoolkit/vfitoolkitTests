function output=DiscP5_downstream(calib)
% P5: feed this block's outputs into an actual toolkit solve
%
% This is the block where the downstream test earns its keep. The two commands return their grids
% in DIFFERENT shapes and in OPPOSITE state orderings (DiscP5_crosstests §1), so a solve is the
% first place either of those can bite. The two variables are given means of opposite sign by the
% calibration - zmean is about [+0.067, -0.067] - so if the toolkit associates the wrong variable
% with the wrong index, the two reported means swap sign and the failure is unmistakable. A
% symmetric calibration would hide exactly this.

fprintf('\n========== P5: downstream (do the outputs go into a solve) ========== \n')

output=struct();
znum=7;
n_d=0; d_grid=[];
n_a=5; a_grid=linspace(0.1,5,n_a)';
n_z=[znum,znum];
Params=struct(); Params.beta=0.9; Params.r=0.03;
DiscountFactorParamNames={'beta'};
ReturnFn=@(aprime,a,z1,z2,r) log(max((1+r)*a+exp(z1)+0.5*exp(z2)-aprime,10^(-10)));
FnsToEvaluate.z1value=@(aprime,a,z1,z2) z1;
FnsToEvaluate.z2value=@(aprime,a,z1,z2) z2;

Mew=calib.full.Mew; Rho=calib.full.Rho; SigmaSq=calib.full.SigmaSq;
zmeanT=calib.full.zmean;
fprintf('the truth is E[z1]=%2.4f and E[z2]=%2.4f - opposite signs, so a swapped ordering is visible \n',zmeanT(1),zmeanT(2))

cmdname={'VAR1_Tauchen (stacked grid)','VAR1_FarmerToda (joint grid)'};
for c_c=1:2
    if c_c==1
        [z_grid,pi_z]=discretizeVAR1_Tauchen(Mew,Rho,SigmaSq,znum,3,struct());
        gridshape='stacked, sum(znum)-by-1';
    else
        opts=struct(); opts.parallel=1; opts.verbose=0;
        [z_grid,pi_z]=discretizeVAR1_FarmerToda(Mew,Rho,SigmaSq,znum,opts);
        gridshape='joint, (znum^M)-by-M';
    end

    vfoptions=struct(); simoptions=struct();
    [V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params,[]);
    AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

    fprintf('%s: grid is %s \n',cmdname{c_c},gridshape)
    fprintf('%s: solve runs, V has the right number of elements [T0], this should be zero: %i \n',cmdname{c_c},numel(V)~=n_a*znum*znum)
    fprintf('%s: StationaryDist sums to one [T1], this should be zero: %2.8e \n',cmdname{c_c},abs(sum(StationaryDist(:))-1))
    fprintf('%s: V is finite everywhere [T0], this should be zero: %i \n',cmdname{c_c},any(~isfinite(V(:))))
    fprintf('%s: AllStats mean of z1 is %+2.4f against %+2.4f [T2], error %2.3e \n',cmdname{c_c},AllStats.z1value.Mean,zmeanT(1),abs(AllStats.z1value.Mean-zmeanT(1)))
    fprintf('%s: AllStats mean of z2 is %+2.4f against %+2.4f [T2], error %2.3e \n',cmdname{c_c},AllStats.z2value.Mean,zmeanT(2),abs(AllStats.z2value.Mean-zmeanT(2)))
    % The explicit swap diagnostic, so a failure says WHY rather than just "the number is wrong"
    swapped=abs(AllStats.z1value.Mean-zmeanT(2))+abs(AllStats.z2value.Mean-zmeanT(1));
    asis=abs(AllStats.z1value.Mean-zmeanT(1))+abs(AllStats.z2value.Mean-zmeanT(2));
    if swapped<asis
        fprintf('%s: THE TWO VARIABLES ARE SWAPPED - the means match the truth better after exchanging them \n',cmdname{c_c})
    end
    output.cmd(c_c).name=cmdname{c_c};
    output.cmd(c_c).means=[AllStats.z1value.Mean,AllStats.z2value.Mean];
    output.cmd(c_c).swapped=(swapped<asis);
end

end
