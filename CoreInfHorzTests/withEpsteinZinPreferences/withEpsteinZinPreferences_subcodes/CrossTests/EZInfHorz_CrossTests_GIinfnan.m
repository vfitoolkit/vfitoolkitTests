function output=EZInfHorz_CrossTests_GIinfnan(n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames)
% DELIBERATE TRIGGER for -Inf in V under the Epstein-Zin grid interpolation layer: the EZ twin of
% CoreInfHorzVFIAlgoTests doPart(11) (InfHorzGI_InfNaN_proposal.md at the toolkit root). The new
% Epstein-Zin postGI raws (2026-10-27) were written WITH the three rules of that proposal built in,
% so there is no red run to keep; this test is what proves the rules are in, and what catches them
% being lost in a later edit.
%
% THE DEFECT CLASS. The postGI Howards step evaluates the chosen fine point as a two-node
% combination lowerprob.*X(lower,:)+(1-lowerprob).*X(lower+1,:). lowerprob is exactly 1 at the
% window's bottom node and exactly 0 at a chosen coarse node, so a zero weight on a node where X is
% -Inf gives 0*(-Inf)=NaN, which the isnan line (meant for pi_z==0) then prices at 0. Under Epstein-Zin
% X is the TRANSFORMED V, (ezc4*V')^ezc5 with -Inf kept as -Inf*ezc4, so the dead node is -Inf or +Inf
% depending on the sign convention, and either one times a zero weight is NaN. The rule: zero each
% weighted term BEFORE the sum where its weight is exactly zero.
%
% WHY THE DEFAULTS CANNOT SEE IT. Howards is switched off after vfoptions.maxhowards (500) iterations
% and the plain iterations that follow converge from the inflated V to the right fixed point. So
% every model is solved twice: at the defaults (census, firing count, V vs ValueFnFromPolicy), and
% with Howards never switched off (maxhowards=maxiter=2000), which on a defective raw cannot
% converge (a correct Bellman step alternating with a wrong Howards step) and is the red/green
% discriminator. Both are compared with pure VFI (howards=0) on the finite set, and the dead set
% with the plain (no-GI) Epstein-Zin solve of the same model.
%
% THE FIRING COUNT, from the policy the solve returns: L2==ngridinterp+2 is a chosen coarse node coded
% as all the weight on the UPPER node, so the zero weight is on the lower node L1; L2==1 is the
% window's bottom node, zero weight on the UPPER node L1+1. A state fires when its zero-weight node
% is -Inf at some z' it can reach. A model where this is 0 is vacuous for this defect and says so.
%
% Two devices, each for all three EZ cases: (1) a_grid extended below the natural borrowing limit
% (dead states BELOW the finite ones: fires the zero weight on a dead LOWER node); (2) a -Inf cap on
% the state at the top z (dead states ABOVE finite ones: fires the zero weight on a dead UPPER node),
% with pi_z altered so the lowest z never reach the top z (else every state is dead).
% Every solve is wrapped in try/catch so a crash prints as a failed check and the bank finishes.

N_z=prod(n_z);
DF=DiscountFactorParamNames;

ReturnFn_cons=@(aprime,a,z,r,w) EZReturnFn_cons_nod_z(aprime,a,z,r,w);
ReturnFn_posU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_positiveUtils_nod_z(aprime,a,z,r,w,ezsigma);
ReturnFn_negU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_negativeUtils_nod_z(aprime,a,z,r,w,ezsigma);
ReturnFn_cons_cap=@(aprime,a,z,r,w,acap,zcap) EZReturnFn_cons_nod_z_cap(aprime,a,z,r,w,acap,zcap);
ReturnFn_posU_cap=@(aprime,a,z,r,w,ezsigma,acap,zcap) EZReturnFn_positiveUtils_nod_z_cap(aprime,a,z,r,w,ezsigma,acap,zcap);
ReturnFn_negU_cap=@(aprime,a,z,r,w,ezsigma,acap,zcap) EZReturnFn_negativeUtils_nod_z_cap(aprime,a,z,r,w,ezsigma,acap,zcap);

vftol=10^(-10);
n2short=5;

fprintf('\n===== Epstein-Zin grid interpolation layer with -Inf in V (deliberate trigger) ===== \n')

%% Device 1: below the natural borrowing limit
n_neg=12;
a_gridNBL=[(-n_neg:1:-1)'; a_grid];
n_aNBL=length(a_gridNBL);
% Device 2: a -Inf cap on the state at the top z
Pcap=Params;
Pcap.acap=0;
Pcap.zcap=(z_grid(end-1)+z_grid(end))/2;
pi_zcap=pi_z;
pi_zcap(1:end-2,end)=0; % the lowest z never reach z_top
pi_zcap=pi_zcap./sum(pi_zcap,2);
fprintf('device 1: a_grid has %i integer-spaced nodes from %i to -1 below the usual cubic grid (n_a=%i); device 2: V(a,z_top)=-Inf for a>%g, pi_z column %i zeroed in rows 1:%i \n',n_neg,-n_neg,n_aNBL,Pcap.acap,N_z,N_z-2)

for device=1:2
    for ezcase=1:3
        vfEZ=struct();
        vfEZ.exoticpreferences='EpsteinZin';
        vfEZ.tolerance=vftol;
        vfEZ.maxiter=2000;
        if ezcase==1
            casestr='cons-units';
            vfEZ.EZutils=0; vfEZ.EZriskaversion='ezgamma'; vfEZ.EZeis='ezphi';
            if device==1, ReturnFn=ReturnFn_cons; else, ReturnFn=ReturnFn_cons_cap; end
        elseif ezcase==2
            casestr='positive utils';
            vfEZ.EZutils=1; vfEZ.EZpositiveutility=1; vfEZ.EZriskaversion='ezrisk';
            if device==1, ReturnFn=ReturnFn_posU; else, ReturnFn=ReturnFn_posU_cap; end
        else
            casestr='negative utils';
            vfEZ.EZutils=1; vfEZ.EZpositiveutility=0; vfEZ.EZriskaversion='ezrisk';
            if device==1, ReturnFn=ReturnFn_negU; else, ReturnFn=ReturnFn_negU_cap; end
        end
        if device==1
            devstr='below NBL'; n_aD=n_aNBL; a_gridD=a_gridNBL; pi_zD=pi_z; PD=Params;
        else
            devstr='cap at top z'; n_aD=n_a; a_gridD=a_grid; pi_zD=pi_zcap; PD=Pcap;
        end
        N_aD=prod(n_aD);
        fprintf('--- %s [EZ %s] --- \n',devstr,casestr)

        % the plain (no-GI) Epstein-Zin solve: the reference dead set
        okp=0; msg='';
        try
            [Vp,~]=ValueFnIter_InfHorz(0,n_aD,n_z,[],a_gridD,z_grid,pi_zD,ReturnFn,PD,DF,[],vfEZ);
            okp=1;
        catch ME
            msg=ME.message;
        end
        fprintf('%s [EZ %s] plain (no GI) solve ran, this should be one: %d \n',devstr,casestr,okp)
        if okp~=1
            fprintf('   error: %s \n',msg)
            continue
        end
        Vp=gather(Vp);
        deadp=(Vp==-Inf);
        fprintf('%s [EZ %s] plain census: finite %i, -Inf %i, +Inf %i, NaN %i of %i \n',devstr,casestr,sum(isfinite(Vp(:))),sum(deadp(:)),sum(Vp(:)==Inf),sum(isnan(Vp(:))),numel(Vp))
        fprintf('%s [EZ %s] plain dead set is non-empty, count should be WELL ABOVE zero: %i \n',devstr,casestr,sum(deadp(:)))

        % GI: pure VFI (howards=0), the reference on the finite set
        vfGI=vfEZ;
        vfGI.gridinterplayer=1;
        vfGI.ngridinterp=n2short;
        vfGI0=vfGI; vfGI0.howards=0;
        ok0=0; msg='';
        try
            [V0,P0]=ValueFnIter_InfHorz(0,n_aD,n_z,[],a_gridD,z_grid,pi_zD,ReturnFn,PD,DF,[],vfGI0);
            ok0=1;
        catch ME
            msg=ME.message;
        end
        fprintf('%s [EZ %s] GI pure-VFI (howards=0) solve ran, this should be one: %d \n',devstr,casestr,ok0)
        if ok0~=1
            fprintf('   error: %s \n',msg)
            continue
        end
        V0=gather(V0);
        fin0=isfinite(V0);
        fprintf('%s [EZ %s] GI pure-VFI: V has no +Inf and no NaN, this should be one: %d \n',devstr,casestr,(sum(V0(:)==Inf)+sum(isnan(V0(:))))==0)
        fprintf('%s [EZ %s] GI pure-VFI: -Inf set equals the plain solve''s, this should be one: %d \n',devstr,casestr,isequal(V0==-Inf,deadp))
        fprintf('%s [EZ %s] GI pure-VFI: finite set is non-empty, this should be one: %d \n',devstr,casestr,sum(fin0(:))>0)

        % GI at the defaults, and with Howards never switched off
        for hw=1:2
            vfH=vfGI;
            if hw==1
                hwstr='Howards at the defaults';
            else
                hwstr='Howards never switched off (maxhowards=maxiter)';
                vfH.maxhowards=vfH.maxiter;
            end
            okh=0; msg='';
            try
                [VH,PH]=ValueFnIter_InfHorz(0,n_aD,n_z,[],a_gridD,z_grid,pi_zD,ReturnFn,PD,DF,[],vfH);
                okh=1;
            catch ME
                msg=ME.message;
            end
            fprintf('%s [EZ %s] GI, %s: solve ran, this should be one: %d \n',devstr,casestr,hwstr,okh)
            if okh~=1
                fprintf('   error: %s \n',msg)
                continue
            end
            VH=gather(VH);
            L1=gather(reshape(PH(1,:,:),[N_aD,N_z]));
            L2=gather(reshape(PH(2,:,:),[N_aD,N_z]));
            fprintf('%s [EZ %s] GI, %s: census: finite %i, -Inf %i, +Inf %i, NaN %i of %i \n',devstr,casestr,hwstr,sum(isfinite(VH(:))),sum(VH(:)==-Inf),sum(VH(:)==Inf),sum(isnan(VH(:))),numel(VH))
            fprintf('%s [EZ %s] GI, %s: V has no +Inf and no NaN, this should be one: %d \n',devstr,casestr,hwstr,(sum(VH(:)==Inf)+sum(isnan(VH(:))))==0)
            fprintf('%s [EZ %s] GI, %s: -Inf set equals the plain solve''s, this should be one: %d \n',devstr,casestr,hwstr,isequal(VH==-Inf,deadp))
            % the firing count: a zero-weight node that is -Inf at a reachable z'
            reach=gather(pi_zD)>0; % [N_z,N_z]: z can reach z'
            zerolower=(L2==n2short+2); % all the weight on the upper node: the lower node L1 has zero weight
            zeroupper=(L2==1);         % the window's bottom node: the upper node L1+1 has zero weight
            fires=false(N_aD,N_z);
            for z_c=1:N_z
                deadlower=any(VH(L1(:,z_c),reach(z_c,:))==-Inf,2);   % L1 node dead at some reachable z'
                deadupper=any(VH(L1(:,z_c)+1,reach(z_c,:))==-Inf,2); % L1+1 node dead at some reachable z'
                fires(:,z_c)=(zerolower(:,z_c)&deadlower)|(zeroupper(:,z_c)&deadupper);
            end
            fprintf('%s [EZ %s] GI, %s: firing count (states whose zero-weight node is -Inf at a reachable z''), should be WELL ABOVE zero: %i \n',devstr,casestr,hwstr,sum(fires(:)))
            fin=fin0&isfinite(VH);
            fprintf('%s [EZ %s] GI, %s: V vs GI pure-VFI on the finite set, this should be zero to convergence tolerance: %.3e \n',devstr,casestr,hwstr,max(abs(VH(fin)-V0(fin))))
            okfp=0; msg='';
            try
                Vfp=gather(ValueFnFromPolicy_InfHorz(PH,0,n_aD,n_z,[],a_gridD,z_grid,pi_zD,ReturnFn,PD,DF,vfH));
                okfp=1;
            catch ME
                msg=ME.message;
            end
            if okfp==1
                finH=isfinite(VH);
                fprintf('%s [EZ %s] GI, %s: ValueFnFromPolicy on the finite set, this should be zero to convergence tolerance: %.3e \n',devstr,casestr,hwstr,max(abs(Vfp(finH)-VH(finH))))
                fprintf('%s [EZ %s] GI, %s: ValueFnFromPolicy -Inf set equals the solver''s, this should be one: %d \n',devstr,casestr,hwstr,isequal(Vfp==-Inf,VH==-Inf))
            else
                fprintf('%s [EZ %s] GI, %s: ValueFnFromPolicy ran, this should be one: 0 \n',devstr,casestr,hwstr)
                fprintf('   error: %s \n',msg)
            end
        end
    end
end

%%
output=struct();

end
