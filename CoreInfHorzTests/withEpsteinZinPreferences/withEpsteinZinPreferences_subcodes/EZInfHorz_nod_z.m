function output=EZInfHorz_nod_z(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,figure_c)
% InfHorz Epstein-Zin, without d, with z, pure discretization (the only InfHorz EZ path).
%
% Raws reached:
%   lowmemory=0 -> ValueFnIter_InfHorz_EpsteinZin_nod_raw
%   lowmemory=1 -> ValueFnIter_InfHorz_EpsteinZin_LowMem_nod_raw
%
% Run three times, once per Epstein-Zin case (cons-units; utility-units positive; utility-units
% negative). The centrepiece for each is the Epstein-Zin recursion recomputed in this file from
% the FORMULA (gamma, phi, risk appear directly; no ezc constants) at the returned Policy, using
% a_grid, z_grid and pi_z without calling the toolkit:
%   cons:  V = [ F^(1-1/phi) + beta*( E_z[V'^(1-gamma)] )^((1-1/phi)/(1-gamma)) ]^(1/(1-1/phi))
%   pos:   V = F + beta*( E_z[V'^(1-risk)] )^(1/(1-risk))
%   neg:   V = F - beta*( E_z[(-V')^(1+risk)] )^(1/(1+risk))
% In infinite horizon this holds to convergence tolerance, not exactly (the solver stops when
% V moves by less than vfoptions.tolerance). Beside it, the SAME identity evaluated at the
% vNM-collapse risk aversion (gamma=1/phi, or risk=0) is printed and must be WELL ABOVE zero:
% that is what makes a wrong exponent, or a dropped sign convention, loud rather than silent.
%
% The vNM reference solve (exoticpreferences not set) is the anti-vacuity yardstick: the EZ
% solution must differ from it, in V and in Policy. For cons-units the comparable reference is
% vNM with period utility x^(1-gamma)/(1-gamma) (the negativeUtils ReturnFn at ezsigma=ezgamma:
% same risk aversion, EIS 1/gamma instead of phi), with its V put into consumption units by
% ((1-gamma)*V)^(1/(1-gamma)). For the utility-units cases it is the same ReturnFn without
% Epstein-Zin.

% Setup vfoptions
n_d=0; d_grid=[]; %#ok<NASGU>

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn_cons=@(aprime,a,z,r,w) EZReturnFn_cons_nod_z(aprime,a,z,r,w);
ReturnFn_posU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_positiveUtils_nod_z(aprime,a,z,r,w,ezsigma);
ReturnFn_negU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_negativeUtils_nod_z(aprime,a,z,r,w,ezsigma);

% Shorthands for the identity checks (plain CPU doubles throughout)
beta=Params.beta;
r=Params.r;
w=Params.w;
ezgamma=Params.ezgamma;
ezphi=Params.ezphi;
ezrisk=Params.ezrisk;
ezsigma=Params.ezsigma;
a_gridcpu=gather(a_grid);
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
    fprintf('\n--- [EZ %s] without d, with z --- \n',casestr)

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
    fprintf('[EZ %s] size(Policy) should be [1 %d %d], is [%s] \n',casestr,N_a,N_z,num2str(size(Policy1)))
    fprintf('[EZ %s] aprime index in range, this should be one: %d \n',casestr,gather(all(Policy1(:)>=1) && all(Policy1(:)<=N_a)))
    fprintf('[EZ %s] V is finite everywhere, this should be one: %d \n',casestr,gather(all(isfinite(V1(:)))))
    if ezcase==3
        fprintf('[EZ %s] V is negative everywhere (negative-valued utility), this should be one: %d \n',casestr,gather(all(V1(:)<0)))
    else
        fprintf('[EZ %s] V is positive everywhere (consumption units / positive-valued utility), this should be one: %d \n',casestr,gather(all(V1(:)>0)))
    end

    %% The Epstein-Zin recursion, recomputed here from the formula at the returned Policy
    Vg=gather(V1);
    aprimeindex=gather(reshape(Policy1,[N_a,N_z]));
    aprimevals=a_gridcpu(aprimeindex);
    c=(1+r)*avals+w*zvals-aprimevals;
    fprintf('[EZ %s] chosen aprime always affordable, this should be one: %d \n',casestr,all(c(:)>0))
    Fpol=-Inf(N_a,N_z);
    if ezcase==1
        Fpol(c>0)=c(c>0);
    elseif ezcase==2
        Fpol(c>0)=((1+c(c>0)).^(1-ezsigma)-1)/(1-ezsigma);
    else
        Fpol(c>0)=(c(c>0).^(1-ezsigma))/(1-ezsigma);
    end
    % E[ transform(V') | z ] at the chosen aprime, and the same at the vNM-collapse risk aversion
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
        RHSalt=(Fpol.^rho+beta*EVpowalt).^(1/rho); % at gamma=1/phi the inner exponent rho/alpha is one
    elseif ezcase==2
        RHS=Fpol+beta*EVpow.^(1/(1-ezrisk));
        RHSalt=Fpol+beta*EVpowalt;
    else
        RHS=Fpol-beta*EVpow.^(1/(1+ezrisk));
        RHSalt=Fpol-beta*EVpowalt;
    end
    resid=Vg-RHS;
    fprintf('[EZ %s] Epstein-Zin recursion recomputed from the formula at Policy, this should be zero to convergence tolerance: %2.8e \n',casestr,max(abs(resid(:))))
    % The same identity at the vNM-collapse risk aversion. If the raws used the wrong power
    % anywhere (or dropped the sign convention of the negative case), THIS is the line that would
    % move towards zero while the one above moved away from it.
    residalt=Vg-RHSalt;
    fprintf('[EZ %s] the SAME identity at the vNM-collapse risk aversion, this should be WELL ABOVE zero: %2.8f \n',casestr,max(abs(residalt(:))))

    %% Optimality: V is the max, so the EZ objective at the vNM reference policy cannot exceed it
    aprimeindexref=gather(reshape(Policyref,[N_a,N_z]));
    cref=(1+r)*avals+w*zvals-a_gridcpu(aprimeindexref);
    Fref=-Inf(N_a,N_z);
    if ezcase==1
        Fref(cref>0)=cref(cref>0);
    elseif ezcase==2
        Fref(cref>0)=((1+cref(cref>0)).^(1-ezsigma)-1)/(1-ezsigma);
    else
        Fref(cref>0)=(cref(cref>0).^(1-ezsigma))/(1-ezsigma);
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
    % This is the regression test for a mis-wired Howards block (the bug class the QH bank found):
    % the Epstein-Zin Howards block re-forms the recursion from Ftemp and is a second,
    % hand-rolled implementation of it.
    vfoptions1_noH=vfoptions;
    vfoptions1_noH.howards=0;
    [V1noH,Policy1noH]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_noH);
    fprintf('[EZ %s] howards=0 (pure VFI), V, this should be zero to convergence tolerance: %.3e \n',casestr,max(abs(V1(:)-V1noH(:))))
    fprintf('[EZ %s] howards=0 (pure VFI), Policy, this should be zero: %.3e \n',casestr,max(abs(Policy1(:)-Policy1noH(:))))

    %% ValueFnFromPolicy_InfHorz: a second, independent reconstruction of V from Policy
    % It re-solves the policy-evaluation fixed point from scratch under the Epstein-Zin
    % recursion, whereas the identity above evaluates the recursion once at the returned V.
    VfromPolicy=ValueFnFromPolicy_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    fprintf('[EZ %s] ValueFnFromPolicy, this should be zero to convergence tolerance: %2.8e \n',casestr,max(abs(VfromPolicy(:)-V1(:))))

    %% Anti-vacuity: Epstein-Zin must matter
    % If the calibration made the EZ and vNM solutions coincide, every check above would pass on
    % a degenerate model. These lines are what make that visible (and are why this bank's z is
    % riskier than the parent bank's).
    fprintf('[EZ %s] Epstein-Zin matters, max|V-Vref| should be WELL ABOVE zero: %2.8f \n',casestr,max(abs(V1(:)-Vrefcomp(:))))
    fprintf('[EZ %s] Epstein-Zin matters, share of states where Policy differs from the vNM reference should be WELL ABOVE zero: %2.4f \n',casestr,gather(mean(Policy1(:)~=Policyref(:))))
    fprintf('[EZ %s] policy is not pinned at the borrowing constraint, share of states with aprime>1 should be WELL ABOVE zero: %2.4f \n',casestr,gather(mean(Policy1(:)>1)))

    %% Figure (one per case: figure_c, 100+figure_c, 200+figure_c)
    fig=figure(figure_c+100*(ezcase-1)); %#ok<NASGU>
    subplot(1,2,1); plot(a_gridcpu,squeeze(gather(PolicyVals1(1,:,zmid))), a_gridcpu,squeeze(gather(PolicyValsref(1,:,zmid))), a_gridcpu,a_gridcpu,'k:')
    title(['[EZ ',casestr,'] aprime policy (mid z)']); xlabel('a'); legend('Epstein-Zin','vNM reference','45 degree','Location','northwest')
    subplot(1,2,2); plot(a_gridcpu,gather(V1(:,zmid)), a_gridcpu,gather(Vrefcomp(:,zmid)))
    title(['[EZ ',casestr,'] value fn (mid z)']); xlabel('a'); legend('Epstein-Zin','vNM reference (same units)','Location','southeast')
end

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
