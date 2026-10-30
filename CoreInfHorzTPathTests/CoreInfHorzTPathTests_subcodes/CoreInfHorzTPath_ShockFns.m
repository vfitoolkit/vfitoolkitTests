function output=CoreInfHorzTPath_ShockFns(T,PricePath,ParamPath,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline)
% Exogenous shocks given by a function (vfoptions/simoptions.ExogShockFn for z, EiidShockFn for e) on an
% InfHorz transition path, each against the same shocks given as explicit grids, which is the oracle:
%   (a) ExogShockFn with constant parameters vs the explicit z_grid
%   (b) ExogShockFn with a parameter (zscale) on the ParamPath vs the explicit [n_z,T] z_grid it implies
%   (c) EiidShockFn with constant parameters vs the explicit e_grid
%   (d) EiidShockFn with a parameter (escale) on the ParamPath vs the explicit [n_e,T] e_grid it implies
% Each is checked on the standalone commands (ValueFnOnTransPath_InfHorz, AgentDistOnTransPath_InfHorz,
% EvalFnOnTransPath_AggVars_InfHorz) and on the general eqm path (shooting, maxiter=2).
% Until 2026-10-27 ExogShockSetup_InfHorz_TPath read an undefined ParamPathStruct in the ParamPath
% branches, so (b) and (d) errored, and (c) returned an empty e grid.
% V_final is the stationary value fn at the final prices and grids, AgentDist_initial the stationary agent
% distribution at the period-1 prices and grids.

n_d=0;
d_grid=[];
n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;

GeneralEqmEqns.CapitalMarket=@(r,assets,earnings) r-0.02-0.01*assets/(1+earnings);
transpathoptions=struct();
transpathoptions.GEnewprice=3; % the shooting algorithm
transpathoptions.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5};
transpathoptions.maxiter=2; % returns the price path after one update, and the general eqm conditions evaluated on it
transpathoptions.verbose=0;
transpathoptions_sa=struct();

% the scale factors along the path: fall from 1.1 to 1 over the first half, so the final grid is the base grid
scalepath=[linspace(1.1,1,floor(T/2)),ones(1,T-floor(T/2))];
Params.zscale=1;
Params.escale=1;
Params.zconst=1.05;
Params.econst=1.05;

for cc=1:4
    vfoptionsF=struct(); simoptionsF=struct(); % the shocks given by a function
    vfoptionsG=struct(); simoptionsG=struct(); % the same shocks as explicit grids
    ParamPathF=ParamPath;
    switch cc
        case 1
            casename='(a) ExogShockFn, constant parameters';
            n_zc=n_z;
            ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
            FnsToEvaluate.assets=@(aprime,a,z) a;
            FnsToEvaluate.earnings=@(aprime,a,z,w) w*z;
            ExogShockFn=@(zconst) deal(z_grid*zconst,pi_z);
            vfoptionsF.ExogShockFn=ExogShockFn; simoptionsF.ExogShockFn=ExogShockFn;
            z_gridG=z_grid*Params.zconst; z_grid1=z_gridG; z_gridT=z_gridG;
            pi_zc=pi_z;
        case 2
            casename='(b) ExogShockFn, parameter on the ParamPath';
            n_zc=n_z;
            ExogShockFn=@(zscale) deal(z_grid*zscale,pi_z);
            vfoptionsF.ExogShockFn=ExogShockFn; simoptionsF.ExogShockFn=ExogShockFn;
            ParamPathF.zscale=scalepath;
            z_gridG=z_grid*scalepath; z_grid1=z_gridG(:,1); z_gridT=z_gridG(:,T); % [n_z,T]
            pi_zc=pi_z;
        case 3
            casename='(c) EiidShockFn, constant parameters';
            n_zc=0; z_gridG=[]; z_grid1=[]; z_gridT=[]; pi_zc=[];
            ReturnFn=@(aprime,a,e,r,w,sigma) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,sigma);
            FnsToEvaluate=struct();
            FnsToEvaluate.assets=@(aprime,a,e) a;
            FnsToEvaluate.earnings=@(aprime,a,e,w) w*e;
            EiidShockFn=@(econst) deal(e_grid*econst,pi_e);
            vfoptionsF.n_e=n_e; vfoptionsF.EiidShockFn=EiidShockFn;
            simoptionsF.n_e=n_e; simoptionsF.EiidShockFn=EiidShockFn;
            e_gridG=e_grid*Params.econst; e_grid1=e_gridG; e_gridT=e_gridG;
        case 4
            casename='(d) EiidShockFn, parameter on the ParamPath';
            EiidShockFn=@(escale) deal(e_grid*escale,pi_e);
            vfoptionsF.n_e=n_e; vfoptionsF.EiidShockFn=EiidShockFn;
            simoptionsF.n_e=n_e; simoptionsF.EiidShockFn=EiidShockFn;
            ParamPathF.escale=scalepath;
            e_gridG=e_grid*scalepath; e_grid1=e_gridG(:,1); e_gridT=e_gridG(:,T); % [n_e,T]
    end
    if cc>=3
        vfoptionsG.n_e=n_e; vfoptionsG.e_grid=e_gridG; vfoptionsG.pi_e=pi_e;
        simoptionsG.n_e=n_e; simoptionsG.e_grid=e_gridG; simoptionsG.pi_e=pi_e;
        vfoptions1=struct('n_e',n_e,'e_grid',e_grid1,'pi_e',pi_e);
        vfoptionsT=struct('n_e',n_e,'e_grid',e_gridT,'pi_e',pi_e);
    else
        vfoptions1=struct();
        vfoptionsT=struct();
    end

    % V_final and AgentDist_initial, from the explicit grids at the two ends of the path
    Params_ii=Params;
    Params_ii.r=PricePath.r(1);
    Params_ii.w=ParamPath.w(1);
    [~,Policy_ii]=ValueFnIter_InfHorz(n_d,n_a,n_zc,d_grid,a_grid,z_grid1,pi_zc,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions1);
    AgentDist_initial=StationaryDist_InfHorz(Policy_ii,n_d,n_a,n_zc,pi_zc,vfoptions1,Params_ii,[]);
    Params_ii.r=PricePath.r(end);
    Params_ii.w=ParamPath.w(end);
    [V_final,Policy_final]=ValueFnIter_InfHorz(n_d,n_a,n_zc,d_grid,a_grid,z_gridT,pi_zc,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptionsT);

    % The standalone commands
    [VPathF,PolicyPathF]=ValueFnOnTransPath_InfHorz(PricePath, ParamPathF, T, V_final, Policy_final, Params, n_d,n_a,n_zc, d_grid,a_grid,z_grid, pi_zc, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptionsF);
    AgentDistPathF=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPathF, PolicyPathF, n_d,n_a,n_zc, pi_zc, T, Params, transpathoptions_sa, simoptionsF);
    AggVarsPathF=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathF, PolicyPathF, PricePath, ParamPathF, Params, T, n_d, n_a, n_zc, d_grid, a_grid,z_grid, transpathoptions_sa, simoptionsF);
    [VPathG,PolicyPathG]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_zc, d_grid,a_grid,z_gridG, pi_zc, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptionsG);
    AgentDistPathG=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPath, PolicyPathG, n_d,n_a,n_zc, pi_zc, T, Params, transpathoptions_sa, simoptionsG);
    AggVarsPathG=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathG, PolicyPathG, PricePath, ParamPath, Params, T, n_d, n_a, n_zc, d_grid, a_grid,z_gridG, transpathoptions_sa, simoptionsG);
    fprintf('Shock fns on the path, %s, vs explicit grids, VPath, this should be zero: %.3e \n',casename,max(abs(gather(VPathF(:))-gather(VPathG(:)))))
    fprintf('Shock fns on the path, %s, vs explicit grids, PolicyPath, this should be zero: %.3e \n',casename,max(abs(gather(PolicyPathF(:))-gather(PolicyPathG(:)))))
    fprintf('Shock fns on the path, %s, vs explicit grids, AgentDistPath, this should be zero: %.3e \n',casename,max(abs(gather(AgentDistPathF(:))-gather(AgentDistPathG(:)))))
    fprintf('Shock fns on the path, %s, vs explicit grids, AggVarsPath earnings, this should be zero: %.3e \n',casename,max(abs(gather(AggVarsPathF.earnings.Mean(:))-gather(AggVarsPathG.earnings.Mean(:)))))
    % The general eqm path
    [PricePathF,GEcondnPathF]=TransitionPath_InfHorz(PricePath, ParamPathF, T, V_final, AgentDist_initial, n_d,n_a,n_zc, d_grid,a_grid,z_grid, pi_zc, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptionsF, vfoptionsF, []);
    [PricePathG,GEcondnPathG]=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, n_d,n_a,n_zc, d_grid,a_grid,z_gridG, pi_zc, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptionsG, vfoptionsG, []);
    GEF=gather(GEcondnPathF.CapitalMarket(:));
    GEG=gather(GEcondnPathG.CapitalMarket(:));
    fprintf('Shock fns on the path, %s, vs explicit grids, general eqm, PricePath r, this should be zero: %.3e \n',casename,max(abs(gather(PricePathF.r(:))-gather(PricePathG.r(:)))))
    fprintf('Shock fns on the path, %s, vs explicit grids, general eqm, GEcondnPath, this should be zero: %.3e \n',casename,max(abs(GEF(1:T-1)-GEG(1:T-1))))
    if cc==2 || cc==4
        % for scale: the explicit grid held at its final value, which is what a shock fn whose
        % parameter on the path was ignored would give
        if cc==2
            [~,GEcondnPathC]=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, n_d,n_a,n_zc, d_grid,a_grid,z_gridT, pi_zc, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, struct(), struct(), []);
        else
            [~,GEcondnPathC]=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, n_d,n_a,n_zc, d_grid,a_grid,z_gridT, pi_zc, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, vfoptionsT, vfoptionsT, []);
        end
        GEC=gather(GEcondnPathC.CapitalMarket(:));
        fprintf('Shock fns on the path, %s, for scale a constant grid would move GEcondnPath by %.3e (must not be zero) \n',casename,max(abs(GEC(1:T-1)-GEG(1:T-1))))
    end
end

output=struct();

end
