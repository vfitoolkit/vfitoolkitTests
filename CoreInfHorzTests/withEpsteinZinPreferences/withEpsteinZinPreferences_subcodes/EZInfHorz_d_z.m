function output=EZInfHorz_d_z(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,figure_c)
% InfHorz Epstein-Zin, with d, with z, pure discretization (the only InfHorz EZ path).
%
% Raws reached:
%   lowmemory=0 -> ValueFnIter_InfHorz_EpsteinZin_raw
%   lowmemory=1 -> ValueFnIter_InfHorz_EpsteinZin_LowMem_raw
%
% There is ONE with-d path. vfoptions.solnmethod='purediscretization_refinement' (the default
% when n_d>0) and 'purediscretization' both reach ValueFnIter_InfHorz_EpsteinZin_raw, which
% takes a flat max over (d,aprime): no Epstein-Zin Refine raw exists. So no 'refine vs plain'
% check is written here; it would be two names for one solve.
%
% Same three cases and same checks as EZInfHorz_nod_z, with the (d,aprime) policy split: the
% Epstein-Zin recursion recomputed from the formula at the returned (d,aprime), the same identity
% at the vNM-collapse risk aversion as the discriminator, optimality against the vNM reference
% policy, lowmemory, howards=0, ValueFnFromPolicy, anti-vacuity.
% The composite good is x=c^varphi*(1-d)^(1-varphi) with c=(1+r)a+w*d*z-aprime, and d<1.

N_d=prod(n_d);
N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn_cons=@(d,aprime,a,z,r,w,varphi) EZReturnFn_cons_d_z(d,aprime,a,z,r,w,varphi);
ReturnFn_posU=@(d,aprime,a,z,r,w,varphi,ezsigma) EZReturnFn_positiveUtils_d_z(d,aprime,a,z,r,w,varphi,ezsigma);
ReturnFn_negU=@(d,aprime,a,z,r,w,varphi,ezsigma) EZReturnFn_negativeUtils_d_z(d,aprime,a,z,r,w,varphi,ezsigma);

% Shorthands for the identity checks (plain CPU doubles throughout)
beta=Params.beta;
r=Params.r;
w=Params.w;
varphi=Params.varphi;
ezgamma=Params.ezgamma;
ezphi=Params.ezphi;
ezrisk=Params.ezrisk;
ezsigma=Params.ezsigma;
a_gridcpu=gather(a_grid);
d_gridcpu=gather(d_grid);
z_gridcpu=gather(z_grid);
pi_zcpu=gather(pi_z);
avals=repmat(a_gridcpu,1,N_z);
zvals=repmat(z_gridcpu',N_a,1);
zmid=ceil(N_z/2);

vftol=10^(-10); % see the driver header on why not the default 1e-9

for ezcase=1:3
    vfoptions=struct();
    vfoptions.exoticpreferences='EpsteinZin';
    vfoptions.tolerance=vftol;
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptions.EZutils=0;
        vfoptions.EZriskaversion='ezgamma';
        vfoptions.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
        Paramsref=Params;
        Paramsref.ezsigma=Params.ezgamma; % vNM with the same risk aversion
        ReturnFnref=ReturnFn_negU;
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptions.EZutils=1;
        vfoptions.EZpositiveutility=1;
        vfoptions.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
        Paramsref=Params;
        ReturnFnref=ReturnFn_posU;
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptions.EZutils=1;
        vfoptions.EZpositiveutility=0;
        vfoptions.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
        Paramsref=Params;
        ReturnFnref=ReturnFn_negU;
    end
    fprintf('\n--- [EZ %s] with d, with z --- \n',casestr)

    %% The vNM reference solve for this case
    vfref=struct(); % deliberately no exoticpreferences field
    vfref.tolerance=vftol;
    [Vref,Policyref]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFnref,Paramsref,DiscountFactorParamNames,[],vfref);
    if ezcase==1
        Vrefcomp=((1-ezgamma)*Vref).^(1/(1-ezgamma)); % into consumption units
    else
        Vrefcomp=Vref;
    end
    PolicyValsref=PolicyInd2Val_InfHorz(Policyref,n_d,n_a,n_z,d_grid,a_grid,vfref);

    %% The Epstein-Zin solve
    [V1,Policy1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    PolicyVals1=PolicyInd2Val_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,vfoptions);
    fprintf('[EZ %s] PolicyInd2Val ran \n',casestr)

    % Shapes and index ranges
    fprintf('[EZ %s] size(V) should be [%d %d], is [%s] \n',casestr,N_a,N_z,num2str(size(V1)))
    fprintf('[EZ %s] size(Policy) should be [2 %d %d], is [%s] \n',casestr,N_a,N_z,num2str(size(Policy1)))
    fprintf('[EZ %s] d index in range, this should be one: %d \n',casestr,gather(all(all(Policy1(1,:)>=1)) && all(all(Policy1(1,:)<=N_d))))
    fprintf('[EZ %s] aprime index in range, this should be one: %d \n',casestr,gather(all(all(Policy1(2,:)>=1)) && all(all(Policy1(2,:)<=N_a))))
    fprintf('[EZ %s] V is finite everywhere, this should be one: %d \n',casestr,gather(all(isfinite(V1(:)))))
    if ezcase==3
        fprintf('[EZ %s] V is negative everywhere (negative-valued utility), this should be one: %d \n',casestr,gather(all(V1(:)<0)))
    else
        fprintf('[EZ %s] V is positive everywhere (consumption units / positive-valued utility), this should be one: %d \n',casestr,gather(all(V1(:)>0)))
    end

    %% The Epstein-Zin recursion, recomputed here from the formula at the returned (d,aprime)
    Vg=gather(V1);
    dindex=gather(reshape(Policy1(1,:,:),[N_a,N_z]));
    aprimeindex=gather(reshape(Policy1(2,:,:),[N_a,N_z]));
    dvals=d_gridcpu(dindex);
    aprimevals=a_gridcpu(aprimeindex);
    c=(1+r)*avals+w*dvals.*zvals-aprimevals;
    ok=(c>0)&(dvals<1);
    fprintf('[EZ %s] chosen (d,aprime) always feasible (c>0, d<1), this should be one: %d \n',casestr,all(ok(:)))
    x=zeros(N_a,N_z);
    x(ok)=(c(ok).^varphi).*((1-dvals(ok)).^(1-varphi));
    Fpol=-Inf(N_a,N_z);
    if ezcase==1
        Fpol(ok)=x(ok);
    elseif ezcase==2
        Fpol(ok)=((1+x(ok)).^(1-ezsigma)-1)/(1-ezsigma);
    else
        Fpol(ok)=(x(ok).^(1-ezsigma))/(1-ezsigma);
    end
    EVpow=zeros(N_a,N_z);
    EVpowalt=zeros(N_a,N_z);
    for z_c=1:N_z
        Vnext=Vg(aprimeindex(:,z_c),:); % N_a-by-N_z: V at the chosen aprime, for each z'
        if ezcase==1
            EVpow(:,z_c)=(Vnext.^(1-ezgamma))*pi_zcpu(z_c,:)';
            EVpowalt(:,z_c)=(Vnext.^(1-1/ezphi))*pi_zcpu(z_c,:)'; % gamma=1/phi, the vNM collapse point
        elseif ezcase==2
            EVpow(:,z_c)=(Vnext.^(1-ezrisk))*pi_zcpu(z_c,:)';
            EVpowalt(:,z_c)=Vnext*pi_zcpu(z_c,:)'; % risk=0, the vNM collapse point
        else
            EVpow(:,z_c)=((-Vnext).^(1+ezrisk))*pi_zcpu(z_c,:)';
            EVpowalt(:,z_c)=(-Vnext)*pi_zcpu(z_c,:)'; % risk=0
        end
    end
    if ezcase==1
        rho=1-1/ezphi;
        alpha=1-ezgamma;
        RHS=(Fpol.^rho+beta*EVpow.^(rho/alpha)).^(1/rho);
        RHSalt=(Fpol.^rho+beta*EVpowalt).^(1/rho);
    elseif ezcase==2
        RHS=Fpol+beta*EVpow.^(1/(1-ezrisk));
        RHSalt=Fpol+beta*EVpowalt;
    else
        RHS=Fpol-beta*EVpow.^(1/(1+ezrisk));
        RHSalt=Fpol-beta*EVpowalt;
    end
    resid=Vg-RHS;
    fprintf('[EZ %s] Epstein-Zin recursion recomputed from the formula at (d,aprime), this should be zero to convergence tolerance: %2.8e \n',casestr,max(abs(resid(:))))
    residalt=Vg-RHSalt;
    fprintf('[EZ %s] the SAME identity at the vNM-collapse risk aversion, this should be WELL ABOVE zero: %2.8f \n',casestr,max(abs(residalt(:))))

    %% Optimality: V is the max, so the EZ objective at the vNM reference policy cannot exceed it
    dindexref=gather(reshape(Policyref(1,:,:),[N_a,N_z]));
    aprimeindexref=gather(reshape(Policyref(2,:,:),[N_a,N_z]));
    dvalsref=d_gridcpu(dindexref);
    cref=(1+r)*avals+w*dvalsref.*zvals-a_gridcpu(aprimeindexref);
    okref=(cref>0)&(dvalsref<1);
    xref=zeros(N_a,N_z);
    xref(okref)=(cref(okref).^varphi).*((1-dvalsref(okref)).^(1-varphi));
    Fref=-Inf(N_a,N_z);
    if ezcase==1
        Fref(okref)=xref(okref);
    elseif ezcase==2
        Fref(okref)=((1+xref(okref)).^(1-ezsigma)-1)/(1-ezsigma);
    else
        Fref(okref)=(xref(okref).^(1-ezsigma))/(1-ezsigma);
    end
    EVpowref=zeros(N_a,N_z);
    for z_c=1:N_z
        Vnext=Vg(aprimeindexref(:,z_c),:);
        if ezcase==1
            EVpowref(:,z_c)=(Vnext.^(1-ezgamma))*pi_zcpu(z_c,:)';
        elseif ezcase==2
            EVpowref(:,z_c)=(Vnext.^(1-ezrisk))*pi_zcpu(z_c,:)';
        else
            EVpowref(:,z_c)=((-Vnext).^(1+ezrisk))*pi_zcpu(z_c,:)';
        end
    end
    if ezcase==1
        RHSref=(Fref.^rho+beta*EVpowref.^(rho/alpha)).^(1/rho);
    elseif ezcase==2
        RHSref=Fref+beta*EVpowref.^(1/(1-ezrisk));
    else
        RHSref=Fref-beta*EVpowref.^(1/(1+ezrisk));
    end
    fprintf('[EZ %s] V >= the EZ objective at the vNM reference policy, this should be one: %d \n',casestr,all(Vg(:)>=RHSref(:)-10^(-8)))

    %% lowmemory=1 builds the return matrix one z at a time; same arithmetic, so same answer
    vfoptions1_lm=vfoptions;
    vfoptions1_lm.lowmemory=1;
    [V1lm,Policy1lm]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_lm);
    fprintf('[EZ %s] lowmemory=1, V, this should be zero: %.3e \n',casestr,max(abs(V1(:)-V1lm(:))))
    fprintf('[EZ %s] lowmemory=1, Policy, this should be zero: %.3e \n',casestr,max(abs(Policy1(:)-Policy1lm(:))))

    %% Howards improvement is only an accelerator, so turning it off must not change anything
    vfoptions1_noH=vfoptions;
    vfoptions1_noH.howards=0;
    [V1noH,Policy1noH]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_noH);
    fprintf('[EZ %s] howards=0 (pure VFI), V, this should be zero to convergence tolerance: %.3e \n',casestr,max(abs(V1(:)-V1noH(:))))
    fprintf('[EZ %s] howards=0 (pure VFI), Policy, this should be zero: %.3e \n',casestr,max(abs(Policy1(:)-Policy1noH(:))))

    %% ValueFnFromPolicy_InfHorz: a second, independent reconstruction of V from Policy
    VfromPolicy=ValueFnFromPolicy_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    fprintf('[EZ %s] ValueFnFromPolicy, this should be zero to convergence tolerance: %2.8e \n',casestr,max(abs(VfromPolicy(:)-V1(:))))

    %% Anti-vacuity: Epstein-Zin must matter
    fprintf('[EZ %s] Epstein-Zin matters, max|V-Vref| should be WELL ABOVE zero: %2.8f \n',casestr,max(abs(V1(:)-Vrefcomp(:))))
    fprintf('[EZ %s] Epstein-Zin matters, share of states where the aprime policy differs from the vNM reference should be WELL ABOVE zero: %2.4f \n',casestr,gather(mean(Policy1(2,:)~=Policyref(2,:))))
    fprintf('[EZ %s] policy is not pinned at the borrowing constraint, share of states with aprime>1 should be WELL ABOVE zero: %2.4f \n',casestr,gather(mean(Policy1(2,:)>1)))
    fprintf('[EZ %s] labor supply is interior, share of states with 1<d index<%d should be WELL ABOVE zero: %2.4f \n',casestr,N_d,gather(mean(Policy1(1,:)>1 & Policy1(1,:)<N_d)))

    %% Figure (one per case: figure_c, 100+figure_c, 200+figure_c)
    fig=figure(figure_c+100*(ezcase-1)); %#ok<NASGU>
    subplot(1,3,1); plot(a_gridcpu,squeeze(gather(PolicyVals1(2,:,zmid))), a_gridcpu,squeeze(gather(PolicyValsref(2,:,zmid))), a_gridcpu,a_gridcpu,'k:')
    title(['[EZ ',casestr,'] aprime policy (mid z)']); xlabel('a'); legend('Epstein-Zin','vNM reference','45 degree','Location','northwest')
    subplot(1,3,2); plot(a_gridcpu,squeeze(gather(PolicyVals1(1,:,zmid))), a_gridcpu,squeeze(gather(PolicyValsref(1,:,zmid))))
    title(['[EZ ',casestr,'] d policy (mid z)']); xlabel('a'); legend('Epstein-Zin','vNM reference','Location','northeast')
    subplot(1,3,3); plot(a_gridcpu,gather(V1(:,zmid)), a_gridcpu,gather(Vrefcomp(:,zmid)))
    title(['[EZ ',casestr,'] value fn (mid z)']); xlabel('a'); legend('Epstein-Zin','vNM reference (same units)','Location','southeast')
end

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
