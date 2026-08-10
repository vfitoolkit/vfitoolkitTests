function output=CoreInfHorzVFIAlgo_ScanPostGIrepeat(n_z,z_grid,pi_z,Params,DiscountFactorParamNames)
% SCAN 2: sweep the postGI accuracy (|preGI - postGI|) over FOUR dimensions on the d, z, noe model:
%   outer:  n_d = 11,21,51,101,201        (d_grid = linspace(0,1,n_d)')
%   middle: n_a = 51,101,201,301,501      (a_grid = 5*linspace(0,1,n_a)'.^3)
%   table:  maxaprimediff = 5,10,20,30,50  (rows)   x   postGIrepeat = 0..5  (cols)
% preGI (full fine-grid search) is the exact reference. ngridinterp=5 throughout.
%
% This shows how the postGI window/repeat requirement depends on BOTH the decision-grid resolution
% (n_d) and the asset-grid resolution (n_a). maxaprimediff is measured in coarse asset cells.
%
% WARNING: this is 5*5*(1+30) = 775 value-function solves, some large (n_d=201, n_a=501). Slow.

DF=DiscountFactorParamNames;
ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);

n_d_list=[11,21,51,101,201];
n_a_list=[51,101,201,301,501];
madlist=[5,10,20,30,50];
repeatlist=[0,1,2,3,4,5];

fprintf('\n================ SCAN 2: maxaprimediff x postGIrepeat, across n_d and n_a (d, z, noe model) ================\n');

for dd=1:length(n_d_list)
    n_d=n_d_list(dd);
    d_grid=linspace(0,1,n_d)';
    fprintf('\n################  n_d = %d  ################\n',n_d);

    for aa=1:length(n_a_list)
        n_a=n_a_list(aa);
        a_grid=5*linspace(0,1,n_a)'.^3;

        % preGI reference (exact full search) for this (n_d,n_a)
        vfo=struct(); vfo.gridinterplayer=1; vfo.ngridinterp=5; vfo.preGI=1;
        [Vpre,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);

        fprintf('\n---- n_d = %d, n_a = %d ----\n',n_d,n_a);
        fprintf('|preGI-postGI|,  rows = maxaprimediff,  cols = postGIrepeat \n');
        fprintf('%18s','');
        for jj=1:length(repeatlist)
            fprintf('%14s',sprintf('repeat=%d',repeatlist(jj)));
        end
        fprintf('\n');
        for ii=1:length(madlist)
            fprintf('  maxaprimediff=%2d:',madlist(ii));
            for jj=1:length(repeatlist)
                if madlist(ii)>floor((n_a-1)/2) % window would exceed the grid (postGI needs n_a>=1+2*maxaprimediff)
                    fprintf('%14s','n/a');
                    continue
                end
                vfo=struct(); vfo.gridinterplayer=1; vfo.ngridinterp=5; vfo.preGI=0;
                vfo.maxaprimediff=madlist(ii); vfo.postGIrepeat=repeatlist(jj);
                [Vpost,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);
                fprintf('%14.8f',max(abs(Vpre(:)-Vpost(:))));
            end
            fprintf('\n');
        end
    end
end

output=struct();

end
