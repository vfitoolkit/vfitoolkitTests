function [LogLij,LogLi,npaths]=EstimFHorzMLE_Enumerate(mu1,AprimeIdx,AprimeProb,pi_z_J,pi_e_J,logEmis,N_a,N_z,N_e,N_j,j0,jlast,conditional)
% The likelihood of each individual by enumerating every latent path through the grid, WITHOUT merging paths that reach the
% same state (merging is exactly what the forward algorithm does, so this is an independent route to the same number).
% Everything is on the cpu and in plain loops.
%
% mu1: N_s-by-1, the distribution at age 1 (N_s=N_a*N_z*N_e, state index a+N_a*(z-1)+N_a*N_z*(e-1))
% AprimeIdx, AprimeProb: N_s-by-N_probs-by-(N_j-1), the next-period a (joint index) and its probability
% pi_z_J: N_z-by-N_z-by-(N_j-1); pi_e_J: N_e-by-N_j (column jj is the distribution of the e realized in period jj)
% logEmis: N_s-by-N-by-N_j, the log emission of each state for each individual at each age (0 where nothing is observed)
% j0, jlast: N-by-1, first and last age of each individual
% conditional: 1 leaves the first age out of LogLi
%
% LogLij(ii,jj)=log L(observations up to jj) - log L(observations up to jj-1), from j0 to jlast (NaN elsewhere)

N=size(logEmis,2);
LogLij=NaN(N,N_j);
LogLi=zeros(N,1);
npaths=0;
N_probs=size(AprimeIdx,2);
for ii=1:N
    s=find(mu1>0);
    w=mu1(s);
    for jj=1:N_j
        if jj>=j0(ii)
            Lbefore=sum(w);
            w=w.*exp(logEmis(s,ii,jj));
            LogLij(ii,jj)=log(sum(w))-log(Lbefore);
        end
        if jj==jlast(ii)
            break
        end
        % every path branches into every (aprime lottery point, zprime, eprime)
        zind=floor(mod(s-1,N_a*N_z)/N_a)+1;
        news=zeros(0,1);
        neww=zeros(0,1);
        for pp=1:N_probs
            ap=AprimeIdx(s,pp,jj);
            pr=AprimeProb(s,pp,jj);
            for zp=1:N_z
                for ep=1:N_e
                    news=[news; ap+N_a*(zp-1)+N_a*N_z*(ep-1)];
                    neww=[neww; w.*pr.*pi_z_J(zind,zp,jj).*pi_e_J(ep,jj+1)];
                end
            end
        end
        keep=(neww>0);
        s=news(keep);
        w=neww(keep);
    end
    npaths=npaths+numel(w);
    temp=LogLij(ii,:);
    if conditional==1
        temp(j0(ii))=NaN;
    end
    LogLi(ii)=sum(temp(~isnan(temp)));
end

end
