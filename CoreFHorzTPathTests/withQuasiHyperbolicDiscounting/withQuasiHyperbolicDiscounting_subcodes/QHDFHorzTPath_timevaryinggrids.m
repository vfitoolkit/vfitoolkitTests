function output=QHDFHorzTPath_timevaryinggrids(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline)
% z and e grids that vary over the transition path, with quasi-hyperbolic discounting (Naive and Sophisticated).
%
% THE METHOD. Two models that are the same model:
%   (A) z_grid as [n_z,T] and e_grid as [n_e,T], drifting up by a factor 1+0.1*t/T
%   (B) the plain z_grid and e_grid, with the same drift entering as parameters zscale and escale on ParamPath:
%       the ReturnFn and the FnsToEvaluate see z*zscale and e*escale
% In both, pi_z and pi_e move over the path (a third of the way towards uniform by period T). So everything
% must agree exactly: (A) is the shock setup's time-varying grids, (B) does not use them at all. This is the
% quasi-hyperbolic value fn path with fastOLG=1, z and e (which read the period's e grid and pi_e with too
% few indices), and the GE path's quasi-hyperbolic value fn step (which read a z grid field that does not exist).
%   (Q1) value fn path, (A) vs (B): V, Policy, Valt and Policyalt, fastOLG=0 and fastOLG=1
%   (Q2) agent dist path, (A) vs (B), fastOLG=0 and fastOLG=1
%   (Q3) one iteration of the GE path (a general eqm condition that depends on earnings), (A) vs (B),
%        fastOLG=0 and fastOLG=1
%   and, for scale, how far the drift moves V and how far one GE iteration moves r (must not be zero)

n_d=0;
d_grid=[];
n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;
driftfactor=1+0.1*(1:1:T)/T; % [1,T]
z_gridT=z_grid.*driftfactor; % [n_z,T]
e_gridT=e_grid.*driftfactor; % [n_e,T]
pimixT=0.3*(1:1:T)/T; % [1,T]
pi_zT=(1-reshape(pimixT,[1,1,T])).*pi_z+reshape(pimixT,[1,1,T]).*ones(n_z,n_z)/n_z; % [n_z,n_z,T]
pi_eT=(1-pimixT).*pi_e+pimixT.*ones(n_e,1)/n_e; % [n_e,T]

ReturnFnA=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluateA=struct();
FnsToEvaluateA.assets=@(aprime,a,z,e) a;
FnsToEvaluateA.earnings=@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;
ReturnFnB=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension,zscale,escale) ReturnFn_nod_z_e_nosemiz(aprime,a,z*zscale,e*escale,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluateB=struct();
FnsToEvaluateB.assets=@(aprime,a,z,e) a;
FnsToEvaluateB.earnings=@(aprime,a,z,e,w,kappa_j,zscale,escale) w*kappa_j*z*zscale*e*escale;
ParamPathB=ParamPath;
ParamPathB.zscale=driftfactor;
ParamPathB.escale=driftfactor;
ParamsB=Params;
ParamsB.zscale=1;
ParamsB.escale=1;

V_final=zeros([n_a,n_z,n_e,N_j],'gpuArray');
Policy_final=ones([1,n_a,n_z,n_e,N_j],'gpuArray');
jequaloneDist=zeros([n_a,n_z,n_e],'gpuArray');
jequaloneDist(1,ceil(n_z/2),ceil(n_e/2))=1;

GeneralEqmEqns=struct();
GeneralEqmEqns.earningsgap=@(earnings) earnings-0.5;

QHtypes={'Naive','Sophisticated'};
for qq=1:2
    vfoptionsA=struct();
    vfoptionsA.exoticpreferences='QuasiHyperbolic';
    vfoptionsA.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
    vfoptionsA.quasi_hyperbolic=QHtypes{qq};
    vfoptionsA.n_e=n_e;
    vfoptionsB=vfoptionsA;
    vfoptionsA.e_grid=e_gridT; vfoptionsA.pi_e=pi_eT;
    vfoptionsB.e_grid=e_grid; vfoptionsB.pi_e=pi_eT;
    simoptionsA=struct();
    simoptionsA.n_e=n_e; simoptionsA.e_grid=e_gridT; simoptionsA.pi_e=pi_eT;
    simoptionsB=struct();
    simoptionsB.n_e=n_e; simoptionsB.e_grid=e_grid; simoptionsB.pi_e=pi_eT;

    % The initial dist (the same for both, from the plain model)
    vfoptions0=vfoptionsA;
    vfoptions0.e_grid=e_grid; vfoptions0.pi_e=pi_e;
    simoptions0=simoptionsA;
    simoptions0.e_grid=e_grid; simoptions0.pi_e=pi_e;
    [~,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFnA,Params,DiscountFactorParamNames,[],vfoptions0);
    AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy0,n_d,n_a,n_z,N_j,pi_z,Params,simoptions0);
    clear Policy0

    %% (Q1), (Q2)
    for ff=0:1 % fastOLG
        transpathoptions=transpathoptionsbaseline;
        transpathoptions.fastOLG=ff;
        simoptionsA.fastOLG=ff;
        simoptionsB.fastOLG=ff;
        % Naive returns (V, Policy, Valt, Policyalt); Sophisticated returns (V, Policy, Vunderbar)
        if qq==1
            [VPathA,PolicyPathA,ValtPathA,PolicyaltPathA]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_gridT, pi_zT, DiscountFactorParamNames, ReturnFnA, transpathoptions, vfoptionsA);
            [VPathB,PolicyPathB,ValtPathB,PolicyaltPathB]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPathB, T, V_final, Policy_final, ParamsB, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_zT, DiscountFactorParamNames, ReturnFnB, transpathoptions, vfoptionsB);
        else
            [VPathA,PolicyPathA,ValtPathA]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_gridT, pi_zT, DiscountFactorParamNames, ReturnFnA, transpathoptions, vfoptionsA);
            [VPathB,PolicyPathB,ValtPathB]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPathB, T, V_final, Policy_final, ParamsB, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_zT, DiscountFactorParamNames, ReturnFnB, transpathoptions, vfoptionsB);
        end
        fprintf('Shocks over the path (QH %s), fastOLG=%i: (Q1) drifting grids vs the drift as parameters, VPath, this should be zero: %.3e \n',QHtypes{qq},ff,gather(max(abs(VPathA(:)-VPathB(:)))))
        fprintf('Shocks over the path (QH %s), fastOLG=%i: (Q1) drifting grids vs the drift as parameters, PolicyPath, this should be zero: %.3e \n',QHtypes{qq},ff,gather(max(abs(PolicyPathA(:)-PolicyPathB(:)))))
        fprintf('Shocks over the path (QH %s), fastOLG=%i: (Q1) drifting grids vs the drift as parameters, ValtPath (Vunderbar for Sophisticated), this should be zero: %.3e \n',QHtypes{qq},ff,gather(max(abs(ValtPathA(:)-ValtPathB(:)))))
        if qq==1
            fprintf('Shocks over the path (QH %s), fastOLG=%i: (Q1) drifting grids vs the drift as parameters, PolicyaltPath, this should be zero: %.3e \n',QHtypes{qq},ff,gather(max(abs(PolicyaltPathA(:)-PolicyaltPathB(:)))))
        end
        AgentDistPathA=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPath, PolicyPathA, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_zT, T,Params, transpathoptions, simoptionsA);
        AgentDistPathB=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPathB, PolicyPathB, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_zT, T,ParamsB, transpathoptions, simoptionsB);
        fprintf('Shocks over the path (QH %s), fastOLG=%i: (Q2) drifting grids vs the drift as parameters, AgentDistPath, this should be zero: %.3e \n',QHtypes{qq},ff,gather(max(abs(AgentDistPathA(:)-AgentDistPathB(:)))))
        if ff==0
            % For scale: the drift must actually move things (the plain model is (B) without the parameters on the path)
            VPathP=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, ParamsB, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_zT, DiscountFactorParamNames, ReturnFnB, transpathoptions, vfoptionsB);
            fprintf('Shocks over the path (QH %s): for scale, the drift moves VPath by up to %.3e (must not be zero, or (Q1)-(Q2) are empty) \n',QHtypes{qq},gather(max(abs(VPathA(:)-VPathP(:)))))
            clear VPathP
        end
        clear VPathA VPathB PolicyPathA PolicyPathB ValtPathA ValtPathB PolicyaltPathA PolicyaltPathB AgentDistPathA AgentDistPathB
    end

    %% (Q3) one iteration of the GE path
    transpathoptionsG=transpathoptionsbaseline;
    transpathoptionsG.maxiter=1;
    transpathoptionsG.GEnewprice=2; % Anderson acceleration
    transpathoptionsG.GEnewprice2.howtoupdate={'earningsgap','r',0,0.1};
    for ff=0:1
        transpathoptionsG.fastOLG=ff;
        simoptionsA.fastOLG=ff;
        simoptionsB.fastOLG=ff;
        PricePathGA=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_gridT, pi_zT, ReturnFnA, FnsToEvaluateA, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsG, simoptionsA, vfoptionsA);
        PricePathGB=TransitionPath_Case1_FHorz(PricePath, ParamPathB, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_zT, ReturnFnB, FnsToEvaluateB, GeneralEqmEqns, ParamsB, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsG, simoptionsB, vfoptionsB);
        fprintf('Shocks over the path (QH %s), fastOLG=%i: (Q3) one iteration of the GE path, drifting grids vs the drift as parameters, this should be zero: %.3e \n',QHtypes{qq},ff,max(abs(PricePathGA.r-PricePathGB.r)))
        fprintf('Shocks over the path (QH %s), fastOLG=%i: for scale, that one iteration moves r by up to %.3e (must not be zero, or (Q3) is empty) \n',QHtypes{qq},ff,max(abs(PricePathGA.r-PricePath.r)))
    end
end

output=struct();

end
