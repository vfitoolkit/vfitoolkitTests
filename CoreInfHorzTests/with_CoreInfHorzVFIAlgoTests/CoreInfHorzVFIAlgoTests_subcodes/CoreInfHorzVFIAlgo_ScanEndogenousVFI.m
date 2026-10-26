function output=CoreInfHorzVFIAlgo_ScanEndogenousVFI(Params,DiscountFactorParamNames)
% SCAN 4 (DIAGNOSTIC): is endogenous VFI (vfoptions.endogenousvfi=1) actually faster?
%
% Endogenous VFI iterates exactly as VFI does, but stops on the span over a of V-Vold within each exogenous state, rather
% than on the sup over all states. That span shrinks at rate beta*alpha_y (alpha_y measures how quickly the endogenous
% state forgets where it started, given the exogenous states), while the sup shrinks only at rate beta. So endogenous VFI
% stops earlier, by how much depending on the model. Nothing else changes, so the only thing to measure is the runtime.
%
% Swept: beta in {0.96, 0.99}, rho_z in {0.9, 0.99}, n_z in {7, 21}, for the z (no e) models with and without d, at
% n_a=500. With Howards at its default, and with Howards off (maxhowards=0, pure VFI), since Howards already does much of
% what endogenous VFI does; and with the grid interpolation layer (ngridinterp=5, Howards at its default). r is set so
% that beta*(1+r)<1 (r=0.8*(1/beta-1)), so assets do not run off the top of the grid.
%
% Nothing here is a check: the V difference is printed against beta/(1-beta)*tolerance, which is the bound on how far each
% of the two solves is from the true V. Speedups are printed rather than asserted.

DF=DiscountFactorParamNames;

beta_list=[0.96,0.99];
rho_list=[0.9,0.99];
n_z_list=[7,21];
n_a=500;
a_grid=5*linspace(0,1,n_a)'.^3;

fprintf('\n================ SCAN 4: endogenous VFI (vfoptions.endogenousvfi=1), runtimes against plain VFI ================\n');

for bb=1:length(beta_list)
    Params.beta=beta_list(bb);
    Params.r=0.8*(1/Params.beta-1);
    for rr=1:length(rho_list)
        rho=rho_list(rr);
        for zz=1:length(n_z_list)
            n_z=n_z_list(zz);
            [z_grid,pi_z]=discretizeAR1_FarmerToda(0,rho,0.1*sqrt(1-rho^2),n_z); % same unconditional std dev of log z for every rho
            z_grid=exp(z_grid);
            for mm=1:2
                if mm==1
                    n_d=0; d_grid=[];
                    ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
                    modelstr='nod, z, noe';
                else
                    n_d=9; d_grid=linspace(0,1,n_d)';
                    ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
                    modelstr='d, z, noe';
                end
                fprintf('\n---- %s, beta=%1.2f, rho_z=%1.2f, n_z=%d, n_a=%d ----\n',modelstr,Params.beta,rho,n_z,n_a);

                vfoptions=struct(); vfoptions.verbose_advice=0;
                Vbound=Params.beta/(1-Params.beta)*10^(-9); % beta/(1-beta)*tolerance, at the default tolerance

                % Warm up (so the first timed solve is not paying GPU/JIT startup costs)
                [~,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions);

                for hh=1:3
                    vfo=vfoptions; vfo.endogenousvfi=0; % plain VFI (the default from beta=0.97 is endogenous VFI)
                    if hh==1
                        hstr='Howards default ';
                    elseif hh==2
                        hstr='Howards OFF     ';
                        vfo.maxhowards=0;
                    else
                        hstr='GI, Howards def ';
                        vfo.gridinterplayer=1; vfo.ngridinterp=5;
                    end
                    tic; [Vp,Policyp]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); wait(gpuDevice); tp=toc;
                    vfo.endogenousvfi=1;
                    tic; [Ve,Policye]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); wait(gpuDevice); te=toc;
                    fprintf('%s plain %2.4fs, endogenous %2.4fs, speedup %2.2f; max|V diff| %.3e (bound per solve %.1e), Pol diff %.3e \n',hstr,tp,te,tp/te,max(abs(Vp(:)-Ve(:))),Vbound,max(abs(Policyp(:)-Policye(:))));
                end
            end
        end
    end
end

%% SCAN 4b: is the gain at high beta just that the Howards defaults were tuned at beta=0.95?
% The vfoptions.howards defaults (40 indexed, 80 sparse) were chosen on this bank at beta=0.95. At beta=0.99 plain VFI may
% simply want more Howards steps, in which case a tuned plain VFI would close the gap to endogenous VFI. So sweep the
% number of Howards steps for both, and compare the BEST plain runtime with the BEST endogenous one. Without GI (indexed
% Howards, the default at these grid sizes) and with GI (sparse Howards, its default). rho_z=0.9, n_z=21, n_a=500.
fprintf('\n================ SCAN 4b: Howards steps x endogenous VFI (does a tuned plain VFI close the gap?) ================\n');
howardslist=[10,20,40,80,160,320];
rho=0.9; n_z=21;
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,rho,0.1*sqrt(1-rho^2),n_z);
z_grid=exp(z_grid);
for bb=1:length(beta_list)
    Params.beta=beta_list(bb);
    Params.r=0.8*(1/Params.beta-1);
    for mm=1:2
        if mm==1
            n_d=0; d_grid=[];
            ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
            modelstr='nod, z, noe';
        else
            n_d=9; d_grid=linspace(0,1,n_d)';
            ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
            modelstr='d, z, noe';
        end
        vfoptions=struct(); vfoptions.verbose_advice=0;
        [~,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions); % warm up
        for gi=0:1
            if gi==0
                gistr='noGI';
            else
                gistr='GI  ';
            end
            fprintf('\n---- %s, beta=%1.2f, rho_z=%1.2f, n_z=%d, n_a=%d, %s ----\n',modelstr,Params.beta,rho,n_z,n_a,gistr);
            tplain=zeros(1,length(howardslist)); tendo=zeros(1,length(howardslist));
            for hh=1:length(howardslist)
                vfo=vfoptions; vfo.howards=howardslist(hh); vfo.endogenousvfi=0; % plain VFI (the default from beta=0.97 is endogenous VFI)
                if gi==1
                    vfo.gridinterplayer=1; vfo.ngridinterp=5;
                end
                tic; [Vp,Policyp]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); wait(gpuDevice); tplain(hh)=toc;
                vfo.endogenousvfi=1;
                tic; [Ve,Policye]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); wait(gpuDevice); tendo(hh)=toc;
                fprintf('%s howards=%3d: plain %2.4fs, endogenous %2.4fs, speedup %2.2f; max|V diff| %.3e, Pol diff %.3e \n',gistr,howardslist(hh),tplain(hh),tendo(hh),tplain(hh)/tendo(hh),max(abs(Vp(:)-Ve(:))),max(abs(Policyp(:)-Policye(:))));
            end
            [tpbest,ip]=min(tplain); [tebest,ie]=min(tendo);
            fprintf('%s BEST: plain %2.4fs (howards=%d), endogenous %2.4fs (howards=%d), best-vs-best speedup %2.2f \n',gistr,tpbest,howardslist(ip),tebest,howardslist(ie),tpbest/tebest);
        end
    end
end

%% SCAN 3d: robustness of the beta-dependent defaults (vfoptions.endogenousvfi from beta>=0.97, and vfoptions.howards)
% The defaults were chosen from SCAN 3c, which used n_a=500 and rho_z=0.9 throughout. Here: a bigger grid (n_a=1500 at
% rho_z=0.9) and a more persistent z (rho_z=0.97 at n_a=500), both at n_z=21, for beta in {0.95,0.96,0.97,0.98,0.99}.
% Each cell times plain and endogenous VFI at every number of Howards steps, and the actual default solve (nothing set),
% and reports the default's runtime against the best of all of them, plus where the threshold would want to be.
fprintf('\n================ SCAN 3d: robustness of the beta-dependent defaults (n_a=1500; rho_z=0.97) ================\n');
howardslist=[20,40,80,160,320,640];
cellslist={[21,1500,0.9],[21,500,0.97]}; % n_z, n_a, rho_z
beta3d_list=[0.95,0.96,0.97,0.98,0.99];
for cc=1:length(cellslist)
    n_z=cellslist{cc}(1); n_a=cellslist{cc}(2); rho=cellslist{cc}(3);
    a_grid=5*linspace(0,1,n_a)'.^3;
    [z_grid,pi_z]=discretizeAR1_FarmerToda(0,rho,0.1*sqrt(1-rho^2),n_z);
    z_grid=exp(z_grid);
    for mm=1:2
        if mm==1
            n_d=0; d_grid=[];
            ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
            modelstr='nod';
        else
            n_d=9; d_grid=linspace(0,1,n_d)';
            ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
            modelstr='d  ';
        end
        for gi=0:1
            vfoptions=struct(); vfoptions.verbose_advice=0;
            if gi==1
                vfoptions.gridinterplayer=1; vfoptions.ngridinterp=5;
                gistr='GI  ';
            else
                gistr='noGI';
            end
            Params.beta=0.95; Params.r=0.8*(1/Params.beta-1);
            [~,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions); % warm up
            for bb=1:length(beta3d_list)
                Params.beta=beta3d_list(bb);
                Params.r=0.8*(1/Params.beta-1);
                th=zeros(1,length(howardslist)); the=zeros(1,length(howardslist));
                for hh=1:length(howardslist)
                    vfo=vfoptions; vfo.howards=howardslist(hh); vfo.endogenousvfi=0;
                    tic; [~,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); wait(gpuDevice); th(hh)=toc;
                    vfo.endogenousvfi=1;
                    tic; [~,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); wait(gpuDevice); the(hh)=toc;
                end
                vfo=vfoptions; % the defaults
                tic; [~,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); wait(gpuDevice); tdef=toc;
                [tbest,ib]=min(th); [tebest,ie]=min(the);
                fprintf('%s %s n_z=%d n_a=%4d rho_z=%1.2f beta=%1.2f: plain %s | endo %s | default %2.4fs = %2.2fx best-of-all | PLAIN best h=%d, ENDO best h=%d, best plain/best endo %2.2f \n',...
                    modelstr,gistr,n_z,n_a,rho,Params.beta,mat2str(th,4),mat2str(the,4),tdef,tdef/min(tbest,tebest),howardslist(ib),howardslist(ie),tbest/tebest);
            end
        end
    end
end

output=struct();

end
