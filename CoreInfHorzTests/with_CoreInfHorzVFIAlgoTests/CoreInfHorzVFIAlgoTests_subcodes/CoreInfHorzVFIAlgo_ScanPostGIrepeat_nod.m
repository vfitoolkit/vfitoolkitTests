function output=CoreInfHorzVFIAlgo_ScanPostGIrepeat_nod(n_z,z_grid,pi_z,Params,DiscountFactorParamNames)
% SCAN 2b: same as CoreInfHorzVFIAlgo_ScanPostGIrepeat, but for the nod, z, noe model.
%   outer: n_a = 51,101,201,301,501       (a_grid = 5*linspace(0,1,n_a)'.^3)
%   table: maxaprimediff = 3,5,8,10,20    (rows)   x   postGIrepeat = 0..5  (cols)
% preGI (full fine-grid search) is the exact reference. ngridinterp=5 throughout.
% There is no n_d loop here (no decision variable), so this is much smaller than the d version.
%
% The point of this scan: the default maxaprimediff for a nod model is 5, and the d-model scan only
% ever tested maxaprimediff values well BELOW its (much wider) defaults. So the maxaprimediff=5 row
% here is the one that says whether postGIrepeat is doing anything at the actual nod default.
% The madlist is therefore centred tightly on 5 rather than spread out as in the d version.
%
% Note: this is 5*(1+30) = 155 value-function solves.

DF=DiscountFactorParamNames;
ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);

n_d=0; d_grid=[];

n_a_list=[51,101,201,301,501];
madlist=[3,5,8,10,20]; % centred on 5, which is the default maxaprimediff for a nod model
repeatlist=[0,1,2,3,4,5];

fprintf('\n================ SCAN 2b: maxaprimediff x postGIrepeat, across n_a (nod, z, noe model) ================\n');
fprintf('[maxaprimediff=5 is the default for a nod model, so that row is the one that matters] \n');

for aa=1:length(n_a_list)
    n_a=n_a_list(aa);
    a_grid=5*linspace(0,1,n_a)'.^3;

    % preGI reference (exact full search) for this n_a
    vfo=struct(); vfo.gridinterplayer=1; vfo.ngridinterp=5; vfo.preGI=1;
    [Vpre,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfo);

    fprintf('\n---- n_a = %d ----\n',n_a);
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

output=struct();

end
