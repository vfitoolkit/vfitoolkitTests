function output=EZFHorzTPath_parampath(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline)
% The Epstein-Zin parameters on the ParamPath. The transition path rebuilds the Epstein-Zin constants
% from each period's parameters, so risk aversion, the EIS and the survival probabilities may all move
% along the path. Two models, nod_z_noe and d_noz_e (so z and e, with and without d), both layouts.
%   (P1)-(P4) a parameter held on the ParamPath at a value DIFFERENT from the one in Params, prices held
%       at their values in Params: every period must equal the stationary solve at the ParamPath value,
%       which (anti-vacuity) must differ from the stationary solve at the Params value. Fails if the path
%       reads the Params value, or reads the ParamPath value once instead of period by period.
%       (P1) EZriskaversion, consumption units; (P2) EZeis, consumption units; (P3) EZriskaversion,
%       negative utils; (P4) the survival probabilities, an age-dependent parameter (N_j-by-T), with
%       warm-glow, consumption units
%   (P5)-(P6) a STEP in the parameter at period s, from the Params value to the new value, with the
%       terminal value fn the stationary solve at the new value. Exactly:
%         every period t>=s equals the stationary solve at the new value
%         every period t<=s-N_j equals the stationary solve at the Params value (everyone alive then is
%           dead before the step reaches them)
%         period s-1 differs from both (should NOT be zero)
%       (P5) EZriskaversion, consumption units; (P6) the survival probabilities, consumption units

PricePathConstant.r=Params.r*ones(1,T);
ParamPathConstant.w=Params.w*ones(1,T);
s=floor(T/2); % the step, for (P5)-(P6)
if s-N_j<1
    error('EZFHorzTPath_parampath needs T>=2*N_j+2 for the step tests')
end

transpathoptionsL=cell(1,2);
transpathoptionsL{1}=transpathoptionsbaseline; transpathoptionsL{1}.fastOLG=0;
transpathoptionsL{2}=transpathoptionsbaseline; transpathoptionsL{2}.fastOLG=1;
layoutname={'slowOLG','fastOLG'};

sj2=[linspace(1,0.8,N_j-1),0]; % different from Params.sj
sjPath2=sj2'.*ones(1,T); % N_j-by-T, an age-dependent parameter on the path

for mm=1:2
    if mm==1
        modelname='nod_z_noe';
        n_dm=0; d_gridm=[]; n_zm=n_z; z_gridm=z_grid; pi_zm=pi_z;
        vfoptionsS=struct();
        ReturnFn_cons=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) EZReturnFn_cons_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,agej,Jr,pension);
        ReturnFn_negU=@(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension) EZReturnFn_negativeUtils_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension);
    else
        modelname='d_noz_e';
        n_dm=n_d; d_gridm=d_grid; n_zm=0; z_gridm=[]; pi_zm=[];
        vfoptionsS=struct();
        vfoptionsS.n_e=vfoptionsbaseline.n_e; vfoptionsS.e_grid=vfoptionsbaseline.e_grid; vfoptionsS.pi_e=vfoptionsbaseline.pi_e;
        ReturnFn_cons=@(d,aprime,a,e,r,w,kappa_j,varphi,agej,Jr,pension) EZReturnFn_cons_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,varphi,agej,Jr,pension);
        ReturnFn_negU=@(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension) EZReturnFn_negativeUtils_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
    end
    fprintf('\n--- EZ-TPath, Epstein-Zin parameters on the ParamPath (%s) ---\n',modelname)
    vfcons=vfoptionsS; vfcons.exoticpreferences='EpsteinZin'; vfcons.EZutils=0; vfcons.EZriskaversion='ezgamma'; vfcons.EZeis='ezphi';
    vfnegU=vfoptionsS; vfnegU.exoticpreferences='EpsteinZin'; vfnegU.EZutils=1; vfnegU.EZpositiveutility=0; vfnegU.EZriskaversion='ezrisk';
    vfconssj=vfcons; vfconssj.survivalprobability='sj'; vfconssj.WarmGlowBequestsFn=@(aprime,wg1,wg2) EZWarmGlowFn_cons(aprime,wg1,wg2);

    %% (P1)-(P4) the parameter held on the ParamPath at a different value
    label={'(P1) EZriskaversion (cons units)','(P2) EZeis (cons units)','(P3) EZriskaversion (negative utils)','(P4) survival probabilities with warm-glow (cons units)'};
    for pp=1:4
        Params2=Params;
        ParamPathP=ParamPathConstant;
        if pp==1
            Params2.ezgamma=4; ParamPathP.ezgamma=4*ones(1,T); vfP=vfcons; ReturnFn=ReturnFn_cons;
        elseif pp==2
            Params2.ezphi=0.7; ParamPathP.ezphi=0.7*ones(1,T); vfP=vfcons; ReturnFn=ReturnFn_cons;
        elseif pp==3
            Params2.ezrisk=5; ParamPathP.ezrisk=5*ones(1,T); vfP=vfnegU; ReturnFn=ReturnFn_negU;
        else
            Params2.sj=sj2; ParamPathP.sj=sjPath2; vfP=vfconssj; ReturnFn=ReturnFn_cons;
        end
        [V_s1,~]=ValueFnIter_Case1_FHorz(n_dm,n_a,n_zm,N_j,d_gridm,a_grid,z_gridm,pi_zm,ReturnFn,Params,DiscountFactorParamNames,[],vfP);
        [V_s2,Policy_s2]=ValueFnIter_Case1_FHorz(n_dm,n_a,n_zm,N_j,d_gridm,a_grid,z_gridm,pi_zm,ReturnFn,Params2,DiscountFactorParamNames,[],vfP);
        fprintf('%s: the stationary solve depends on the parameter, should NOT be zero: %.3e \n',label{pp},max(abs(V_s1(:)-V_s2(:))))
        Vrep=reshape(repmat(V_s2,[ones(1,ndims(V_s2)),T]),[],1);
        Policyrep=reshape(repmat(Policy_s2,[ones(1,ndims(Policy_s2)),T]),[],1);
        for ll=1:2
            [VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePathConstant, ParamPathP, T, V_s2, Policy_s2, Params, n_dm, n_a, n_zm, N_j, d_gridm, a_grid,z_gridm, pi_zm, DiscountFactorParamNames, ReturnFn, transpathoptionsL{ll}, vfP);
            fprintf('%s on the ParamPath, every period against the stationary solve at the ParamPath value (%s), this should be zero, V: %.3e \n',label{pp},layoutname{ll},max(abs(VPath(:)-Vrep)))
            fprintf('%s on the ParamPath, every period against the stationary solve at the ParamPath value (%s), this should be zero, Policy: %.3e \n',label{pp},layoutname{ll},max(abs(PolicyPath(:)-Policyrep)))
        end
        clear VPath PolicyPath Vrep Policyrep V_s1 V_s2 Policy_s2
    end

    %% (P5)-(P6) a step in the parameter at period s
    labelS={'(P5) step in EZriskaversion (cons units)','(P6) step in the survival probabilities (cons units, with warm-glow)'};
    for pp=1:2
        Params2=Params;
        ParamPathS=ParamPathConstant;
        if pp==1
            Params2.ezgamma=4; ParamPathS.ezgamma=[Params.ezgamma*ones(1,s-1),4*ones(1,T-s+1)]; vfP=vfcons;
        else
            Params2.sj=sj2; ParamPathS.sj=[Params.sj'.*ones(1,s-1),sj2'.*ones(1,T-s+1)]; vfP=vfconssj;
        end
        [V_s1,Policy_s1]=ValueFnIter_Case1_FHorz(n_dm,n_a,n_zm,N_j,d_gridm,a_grid,z_gridm,pi_zm,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfP);
        [V_s2,Policy_s2]=ValueFnIter_Case1_FHorz(n_dm,n_a,n_zm,N_j,d_gridm,a_grid,z_gridm,pi_zm,ReturnFn_cons,Params2,DiscountFactorParamNames,[],vfP);
        for ll=1:2
            [VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePathConstant, ParamPathS, T, V_s2, Policy_s2, Params, n_dm, n_a, n_zm, N_j, d_gridm, a_grid,z_gridm, pi_zm, DiscountFactorParamNames, ReturnFn_cons, transpathoptionsL{ll}, vfP);
            VPath=reshape(VPath,[],T);
            PolicyPath=reshape(PolicyPath,[],T);
            fprintf('%s, periods from s on against the stationary solve at the new value (%s), this should be zero, V: %.3e \n',labelS{pp},layoutname{ll},max(max(abs(VPath(:,s:T)-V_s2(:)))))
            fprintf('%s, periods from s on against the stationary solve at the new value (%s), this should be zero, Policy: %.3e \n',labelS{pp},layoutname{ll},max(max(abs(PolicyPath(:,s:T)-Policy_s2(:)))))
            fprintf('%s, periods up to s-N_j against the stationary solve at the Params value (%s), this should be zero, V: %.3e \n',labelS{pp},layoutname{ll},max(max(abs(VPath(:,1:s-N_j)-V_s1(:)))))
            fprintf('%s, periods up to s-N_j against the stationary solve at the Params value (%s), this should be zero, Policy: %.3e \n',labelS{pp},layoutname{ll},max(max(abs(PolicyPath(:,1:s-N_j)-Policy_s1(:)))))
            fprintf('%s, period s-1 against the stationary solve at the new value (%s), should NOT be zero: %.3e \n',labelS{pp},layoutname{ll},max(abs(VPath(:,s-1)-V_s2(:))))
            fprintf('%s, period s-1 against the stationary solve at the Params value (%s), should NOT be zero: %.3e \n',labelS{pp},layoutname{ll},max(abs(VPath(:,s-1)-V_s1(:))))
        end
        clear VPath PolicyPath V_s1 V_s2 Policy_s1 Policy_s2
    end
end

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
