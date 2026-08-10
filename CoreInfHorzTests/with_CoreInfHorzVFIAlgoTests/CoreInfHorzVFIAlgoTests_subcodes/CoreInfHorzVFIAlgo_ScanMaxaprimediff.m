function output=CoreInfHorzVFIAlgo_ScanMaxaprimediff(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames)
% SCAN 1: on the d, z, noe model, vary postGI maxaprimediff = 5,8,10,15,20 and report
% |preGI - postGI| for each. preGI (full fine-grid search) is the exact reference, so this
% shows how wide the postGI search window must be before it stops clipping the optimum
% (i.e. where |preGI-postGI| hits ~0). ngridinterp=5 throughout.

DF=DiscountFactorParamNames;
ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);

fprintf('\n================ SCAN 1: postGI maxaprimediff (d, z, noe model) ================\n');

% preGI reference (exact full search)
vfo=struct(); vfo.gridinterplayer=1; vfo.ngridinterp=5; vfo.preGI=1;
[Vpre,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);

madlist=[5,8,10,15,20];
for ii=1:length(madlist)
    vfo=struct(); vfo.gridinterplayer=1; vfo.ngridinterp=5; vfo.preGI=0; vfo.maxaprimediff=madlist(ii);
    [Vpost,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
    fprintf('maxaprimediff=%2d: |preGI-postGI| = %2.8f \n',madlist(ii),max(abs(Vpre(:)-Vpost(:))));
end

output=struct();

end
