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
% what endogenous VFI does. r is set so that beta*(1+r)<1 (r=0.8*(1/beta-1)), so assets do not run off the top of the grid.
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

                for hh=1:2
                    vfo=vfoptions;
                    if hh==1
                        hstr='Howards default ';
                    else
                        hstr='Howards OFF     ';
                        vfo.maxhowards=0;
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

output=struct();

end
