function output=DiscP9_MVNormal(calib,figure_c)
% P9: MVNormal_ProbabilitiesOnGrid
%
% Before this subcode the command was called from exactly two places in the bank, both at M=1, both
% guarding the B12 regression. Everything multivariate about it - which is all of it, given the name
% - was unexercised. So the invariants come first and the accuracy second.
%
% THE OUTPUT SHAPE IS NOT UNIFORM ACROSS M, and that is the first thing to pin. At l_z==1 the command
% returns a column vector of length znum; at l_z>=2 it reshapes to an array of size znum, so a
% 7-by-7 at M=2 and a 2-by-2-by-2-by-2-by-2 at M=5. A caller that wrote size(P,1) would be right at
% M=1 and wrong everywhere else. Both shapes are asserted below rather than assumed.
%
% AND THE PROBABILITIES ARE NOT GUARANTEED TO SUM TO ONE. The bins tile the whole of R^M - the
% outermost spacings are Inf - so mathematically they must. Numerically the l_z>=2 branch goes
% through mvncdf, which is accurate only to a tolerance, and B7's implementation notes that it
% renormalises for exactly this reason. So the sum is measured and reported against a loose bar, not
% asserted as an identity, and how far it drifts with M is itself the finding.

fprintf('\n========== P9: MVNormal_ProbabilitiesOnGrid ========== \n')

output=struct();
nc=length(calib.mv);
sumerr=zeros(1,nc); meanerr=zeros(1,nc); varerr=zeros(1,nc); coverr=zeros(1,nc);

for c_c=1:nc
    M=calib.mv(c_c).M; znum=calib.mv(c_c).znum;
    Mew=calib.mv(c_c).Mew; Sigma=calib.mv(c_c).Sigma; z_grid=calib.mv(c_c).z_grid;
    mvo=struct(); mvo.parallel=1; mvo.verbose=0;

    P=MVNormal_ProbabilitiesOnGrid(z_grid,Mew,Sigma,znum,mvo);
    P=gather(P);

    % --- shape, which differs between M=1 and M>1
    if M==1
        fprintf('M=1: P is a column vector of length znum [T0], this should be zero: %i \n',any(size(P)~=[znum(1),1]))
    else
        fprintf('M=%i: P is an array of size znum [T0], this should be zero: %i \n',M,any(size(P)~=znum))
    end
    fprintf('M=%i: P has prod(znum) elements [T0], this should be zero: %i \n',M,numel(P)~=prod(znum))
    fprintf('M=%i: P is non-negative [T0], this should be zero: %i \n',M,any(P(:)<0))
    fprintf('M=%i: P has no NaN or Inf [T0], this should be zero: %i \n',M,any(~isfinite(P(:))))
    sumerr(c_c)=abs(sum(P(:))-1);
    % The bar is set from the 2026-09-16 run, the first time this was measured. The drift was
    % exactly 0 at M=1 and M=2, 6.7e-16 at M=3, then 8.8e-06 at M=4 and 9.8e-05 at M=5 - flat at
    % machine precision until it is not, and then two orders per added dimension. That is mvncdf's
    % tolerance compounding over more rectangles, and it is why the callers of this command
    % renormalise (see B7). A single bar across M would be vacuous at the bottom and failing at the
    % top, so there is one per M, each about an order above what was measured. The grid sizes moved
    % at M=4 and M=5 after that run, so these two may need resetting once they have been re-measured.
    sumbar=[1e-10,1e-10,1e-10,1e-04,1e-03];
    fprintf('M=%i: P sums to one, which is exact in the maths and only approximate through mvncdf [T2], this should be below %g: %2.3e \n',M,sumbar(c_c),sumerr(c_c))

    % --- moments, against the distribution it is approximating
    % Renormalise first, which is what the callers of this command do (see B7). Without it the moment
    % errors would be contaminated by the sum drift measured just above, and the two would not be
    % separable.
    Pv=P(:)/sum(P(:));
    zvals=CreateGridvals(znum(:),z_grid,1);
    mhat=(Pv'*zvals)';
    dev=zvals-mhat';
    Vhat=(dev.*Pv)'*dev;
    meanerr(c_c)=max(abs(mhat-Mew));
    varerr(c_c)=max(abs(diag(Vhat)-diag(Sigma)));
    offd=~logical(eye(M));
    if M>1
        coverr(c_c)=max(abs(Vhat(offd)-Sigma(offd)));
    else
        coverr(c_c)=0;
    end
    fprintf('M=%i: worst mean error [T2] %2.3e, worst variance error [T2] %2.3e, worst covariance error [T2] %2.3e \n',M,meanerr(c_c),varerr(c_c),coverr(c_c))

    % --- the off-diagonals are actually being used
    % With the same variances but zero covariances the answer must be DIFFERENT, or the command is
    % ignoring the off-diagonal entries of Sigma - which at M=1 is vacuous and at M>1 is the whole
    % point of using a multivariate command.
    if M>1
        Pd=gather(MVNormal_ProbabilitiesOnGrid(z_grid,Mew,calib.mv(c_c).Sigmadiag,znum,mvo));
        fprintf('M=%i: zeroing the off-diagonals of Sigma changes P, so they are read [T0], this should be one: %i \n',M,max(abs(Pd(:)-P(:)))>1e-08)
        % ...and with zero covariances the joint must factor into its marginals.
        Pdv=Pd(:)/sum(Pd(:));
        Pdarr=reshape(Pdv,[znum(:)',1]);
        marg=cell(1,M);
        for m_c=1:M
            A=Pdarr;
            for k_c=M:-1:1
                if k_c~=m_c
                    A=sum(A,k_c);
                end
            end
            marg{m_c}=A(:);
        end
        Pfac=marg{1};
        for m_c=2:M
            Pfac=kron(marg{m_c},Pfac); % variable 1 fastest, matching CreateGridvals
        end
        fprintf('M=%i: with zero covariances the joint is the product of its marginals [T1], this should be zero: %2.8e \n',M,max(abs(Pfac-Pdv)))
    end

    % --- Sigma as a vector is the documented shorthand for a diagonal matrix
    Pvec=gather(MVNormal_ProbabilitiesOnGrid(z_grid,Mew,diag(calib.mv(c_c).Sigmadiag),znum,mvo));
    Pmat=gather(MVNormal_ProbabilitiesOnGrid(z_grid,Mew,calib.mv(c_c).Sigmadiag,znum,mvo));
    fprintf('M=%i: Sigma given as a vector equals the same Sigma given as a diagonal matrix [T1], this should be zero: %2.8e \n',M,max(abs(Pvec(:)-Pmat(:))))
end

fprintf('\nthe sum-to-one drift across M=1 to %i:',nc);
for c_c=1:nc
    fprintf(' %2.1e',sumerr(c_c));
end
fprintf(' \n   (M=1 goes through erfc and should be at machine precision; M>=2 goes through mvncdf \n')
fprintf('   and should not. If the drift grows with M that is mvncdf''s tolerance compounding over \n')
fprintf('   more dimensions, which is a reason to renormalise rather than a defect.) \n')

%% Error paths
mvo=struct(); mvo.parallel=1; mvo.verbose=0;
% B13 widened the guard from l_z>=5 to l_z>5, which is what made the five-dimension case above
% reachable. Six must still be rejected, or the guard went one too far.
M6=6; zn6=2*ones(1,M6); zg6=repmat(linspace(-1,1,2)',M6,1);
try
    MVNormal_ProbabilitiesOnGrid(zg6,zeros(M6,1),eye(M6),zn6,mvo);
    fprintf('\nsix dimensions was accepted, but the command documents a maximum of five [T0], this should be one: 0 \n')
catch ME
    fprintf('\nsix dimensions errors, and the message says five, this should be one: %i \n',contains(ME.message,'five'))
end
try
    MVNormal_ProbabilitiesOnGrid([1;2;3],0,1,5,mvo);
    fprintf('a z_grid of the wrong length was accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('a z_grid of the wrong length errors, and the message names z_grid, this should be one: %i \n',contains(ME.message,'z_grid'))
end

output.sumerr=sumerr; output.meanerr=meanerr; output.varerr=varerr; output.coverr=coverr;

%% Figure
figure(figure_c)
Ms=1:nc;
subplot(1,2,1)
semilogy(Ms,max(sumerr,1e-18),'o-',Ms,max(meanerr,1e-18),'s-',Ms,max(varerr,1e-18),'d-',Ms,max(coverr,1e-18),'^-')
xlabel('M (number of dimensions)'); ylabel('|error|')
title('MVNormal\_ProbabilitiesOnGrid'); legend('sum to one','mean','variance','covariance','Location','best')
subplot(1,2,2)
P1=gather(MVNormal_ProbabilitiesOnGrid(calib.mv(1).z_grid,calib.mv(1).Mew,calib.mv(1).Sigma,calib.mv(1).znum,mvo));
plot(calib.mv(1).z_grid,P1,'o-')
xlabel('z'); ylabel('probability'); title('M=1: probabilities on the given grid')
sgtitle('P9: MVNormal\_ProbabilitiesOnGrid')

end
