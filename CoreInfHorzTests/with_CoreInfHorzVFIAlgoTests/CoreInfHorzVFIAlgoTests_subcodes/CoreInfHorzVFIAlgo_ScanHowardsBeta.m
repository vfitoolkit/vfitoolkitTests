function output=CoreInfHorzVFIAlgo_ScanHowardsBeta(Params,DiscountFactorParamNames)
% SCAN 3c (DIAGNOSTIC): how should the default number of Howards steps (vfoptions.howards) depend on beta, and does a
% tuned plain VFI close the gap to endogenous VFI (vfoptions.endogenousvfi=1) at high beta?
%
% The defaults (40 with indexed Howards, 80 with sparse Howards) were chosen on this bank at beta=0.95. The endogenous VFI
% scan (SCAN 4b) found that at beta=0.99 plain VFI is fastest at 160-320 Howards steps, and that the defaults cost it up
% to about 1.5x. The level of V converges at rate beta, and each Howards step is cheap relative to a Bellman step, so the
% right number of Howards steps plausibly grows like 1/(1-beta). This scan measures it.
%
% Each cell is solved by plain VFI and by endogenous VFI at every number of Howards steps; the line reports the best of
% each, the cost of the current default, and best-plain against best-endogenous (which says from which beta, if any,
% endogenous VFI is worth turning on by default). Swept: beta in {0.9, 0.95, 0.97, 0.98, 0.99, 0.995}, vfoptions.howards in
% {20, 40, 80, 160, 320, 640}, with and without d, with and without the grid interpolation layer, and (n_z,n_a) in
% {(7,500), (21,500), (75,500)}. The Howards method is left at its default, which is indexed at n_z=7 and 21 without GI
% and sparse at n_z=75 or with GI, so both methods are covered. r=0.8*(1/beta-1) so that beta*(1+r)<1. rho_z=0.9.
%
% One line per cell gives every runtime, the best howards, and the cost of the current default relative to the best.
% Nothing here is a check.

DF=DiscountFactorParamNames;

beta_list=[0.9,0.95,0.97,0.98,0.99,0.995];
howardslist=[20,40,80,160,320,640];
nz_list=[7,21,75];
n_a=500;
a_grid=5*linspace(0,1,n_a)'.^3;
rho=0.9;

fprintf('\n================ SCAN 3c: number of Howards steps against beta (plain VFI) ================\n');
fprintf('howards list: %s \n',mat2str(howardslist));

for zz=1:length(nz_list)
    n_z=nz_list(zz);
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
            % the Howards method the default would pick (see ValueFnIter_InfHorz)
            if gi==1 || n_z>=50
                sparsestr='sparse ';
                howardsdefault=80;
            else
                sparsestr='indexed';
                howardsdefault=40;
            end
            Params.beta=0.95; Params.r=0.8*(1/Params.beta-1);
            [~,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions); % warm up
            for bb=1:length(beta_list)
                Params.beta=beta_list(bb);
                Params.r=0.8*(1/Params.beta-1);
                th=zeros(1,length(howardslist)); the=zeros(1,length(howardslist));
                for hh=1:length(howardslist)
                    vfo=vfoptions; vfo.howards=howardslist(hh);
                    tic; [~,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); wait(gpuDevice); th(hh)=toc;
                    vfo.endogenousvfi=1;
                    tic; [~,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo); wait(gpuDevice); the(hh)=toc;
                end
                [tbest,ib]=min(th);
                [tebest,ie]=min(the);
                tdefault=th(howardslist==howardsdefault);
                fprintf('%s %s n_z=%3d %s beta=%1.3f: plain %s | endo %s | PLAIN best howards=%3d (%2.4fs), default %d costs %2.2fx | ENDO best howards=%3d (%2.4fs), endo default costs %2.2fx | best plain/best endo %2.2f \n',...
                    modelstr,gistr,n_z,sparsestr,Params.beta,mat2str(th,4),mat2str(the,4),howardslist(ib),tbest,howardsdefault,tdefault/tbest,howardslist(ie),tebest,the(howardslist==howardsdefault)/tebest,tbest/tebest);
            end
        end
    end
end

output=struct();

end
