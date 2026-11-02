function output=EZFHorzTPath_refusals(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline)
% What Epstein-Zin on the FHorz transition path does NOT do must be refused by name, at the two top-level
% commands (ValueFnOnTransPath_Case1_FHorz and TransitionPath_Case1_FHorz), never silently solved as
% something else. Checked AS A SET: every line is 'should be one', so a refusal that stops firing, or one
% that fires with a different message, is a failure rather than a vanished check. Each refusal must name
% the cause; the error message is printed alongside.
%   (R1) a semi-exogenous state
%   (R2) an experience asset
%   (R3) two endogenous states
%   (R4) a missing vfoptions.EZriskaversion; (R5) EZutils=0 without vfoptions.EZeis
%   (R6) a survival-probability parameter that is not of length N_j

DF=DiscountFactorParamNames;
ReturnFn=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) EZReturnFn_cons_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,agej,Jr,pension);
FnsToEvaluate.assets=@(aprime,a,z) a;
GeneralEqmEqns.CapitalMarket=@(r,assets) r-(0.3*((assets+0.1)^(-0.7))-0.08);

vfEZ=struct();
vfEZ.exoticpreferences='EpsteinZin';
vfEZ.EZutils=0; vfEZ.EZriskaversion='ezgamma'; vfEZ.EZeis='ezphi';

% Terminal objects and an initial agent distribution of the right shapes. The refusals must fire before
% any of them is used.
V_final=zeros([n_a,n_z,N_j],'gpuArray');
Policy_final=ones([1,n_a,n_z,N_j],'gpuArray');
jequaloneDist=zeros([n_a,n_z],'gpuArray'); jequaloneDist(1,ceil(n_z/2))=1;
AgentDist_initial=ones([n_a,n_z,N_j],'gpuArray')/(n_a*n_z*N_j);

transpathoptionsGE=transpathoptionsbaseline;
transpathoptionsGE.GEnewprice=3;
transpathoptionsGE.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1};
transpathoptionsGE.maxiter=1;
transpathoptionsGE.verbose=0;

refusalname={'(R1) semi-exogenous state','(R2) experience asset','(R3) two endogenous states','(R4) missing EZriskaversion','(R5) EZutils=0 without EZeis','(R6) survival probabilities not of length N_j'};
refusalword={'semi-exogenous','experience','endogenous states','EZriskaversion','EZeis','N_j'};
for rr=1:6
    vfR=vfEZ;
    n_aR=n_a; a_gridR=a_grid; ParamsR=Params;
    V_finalR=V_final; Policy_finalR=Policy_final; jequaloneDistR=jequaloneDist; AgentDist_initialR=AgentDist_initial;
    if rr==1
        vfR.n_semiz=vfoptionsbaseline.n_semiz; vfR.semiz_grid=vfoptionsbaseline.semiz_grid; vfR.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
    elseif rr==2
        vfR.experienceasset=1;
        vfR.aprimeFn=@(aprime,a2) a2;
    elseif rr==3
        n_aR=[n_a,2]; a_gridR=[a_grid; 0; 1];
        V_finalR=zeros([n_a,2,n_z,N_j],'gpuArray'); Policy_finalR=ones([2,n_a,2,n_z,N_j],'gpuArray');
        jequaloneDistR=zeros([n_a,2,n_z],'gpuArray'); jequaloneDistR(1,1,ceil(n_z/2))=1;
        AgentDist_initialR=ones([n_a,2,n_z,N_j],'gpuArray')/(n_a*2*n_z*N_j);
    elseif rr==4
        vfR=rmfield(vfR,'EZriskaversion');
    elseif rr==5
        vfR=rmfield(vfR,'EZeis');
    else
        vfR.survivalprobability='sjbad';
        ParamsR.sjbad=ones(1,N_j-1);
    end
    % ValueFnOnTransPath_Case1_FHorz, both layouts
    for fastOLG=0:1
        transpathoptionsR=transpathoptionsbaseline; transpathoptionsR.fastOLG=fastOLG;
        fired=0; msg='(it ran)';
        try
            ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_finalR, Policy_finalR, ParamsR, 0, n_aR, n_z, N_j, [], a_gridR,z_grid, pi_z, DF, ReturnFn, transpathoptionsR, vfR);
        catch ME
            msg=ME.message;
            fired=contains(msg,refusalword{rr});
        end
        fprintf('%s with Epstein-Zin: ValueFnOnTransPath_Case1_FHorz (fastOLG=%i) refuses it by name, this should be one: %i \n',refusalname{rr},fastOLG,fired)
        fprintf('   message: %s \n',msg)
    end
    % TransitionPath_Case1_FHorz
    fired=0; msg='(it ran)';
    try
        TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_finalR, AgentDist_initialR, jequaloneDistR, 0, n_aR, n_z, N_j, [], a_gridR,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsR, DF, AgeWeightParamNames, transpathoptionsGE, struct(), vfR);
    catch ME
        msg=ME.message;
        fired=contains(msg,refusalword{rr});
    end
    fprintf('%s with Epstein-Zin: TransitionPath_Case1_FHorz refuses it by name, this should be one: %i \n',refusalname{rr},fired)
    fprintf('   message: %s \n',msg)
end

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
