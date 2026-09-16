% Setup for P9: a normal distribution on a GIVEN grid
%
% P9 is not a process class. There is no dynamics here - no transition matrix, no persistence, no
% ages - just the question of how to allocate probabilities to grid points you were handed rather
% than got to choose. It sits at the end of the process blocks because that is where it belongs in
% the reading order, not because it is one of them.
%
% WHAT MAKES IT A BLOCK RATHER THAN A ONE-COMMAND FORMALITY. Three commands can now be asked to put
% a normal distribution on the same grid:
%    MVNormal_ProbabilitiesOnGrid                    natively, at any M from 1 to 5
%    discretizeIIDNormal_Tauchen    via e_grid       added for this block
%    discretizeIIDNormal_TanakaToda via e_grid       added for this block
% At M=1 all three are answering exactly the same question on exactly the same grid, which turns
% "probabilities on a grid" into a genuine method comparison - and the three do NOT agree, because
% they are not trying to do the same thing. Tauchen and MVNormal both integrate the density over the
% bin around each point, so they should agree to machine precision. Tanaka-Toda instead solves a
% maximum entropy problem that MATCHES MOMENTS on that grid, so it is deliberately different, and on
% a grid wide enough to carry them it should be MORE accurate on the moments and less accurate as an
% approximation to the bin probabilities. Saying which is "right" depends on what the grid is for.
%
% TWO THINGS HERE HAVE NEVER RUN. MVNormal_ProbabilitiesOnGrid is called from exactly two places in
% the bank before this block, both at M=1, guarding the B12 regression (it used to pass a variance
% where the normal cdf wanted a standard deviation). Everything multivariate about it is unexercised.
% And B13 widened its dimension guard from l_z>=5 to l_z>5, which made the l_z==5 branch reachable
% for the first time - the proposal's note is that it "has presumably never run". This block is
% where it runs.

calibNG=struct();

%% The distribution
calibNG.mew=0.3;          % deliberately not zero, so a command that ignored the mean would show
calibNG.sigma=0.7;
calibNG.var=calibNG.sigma^2;

%% The grids, at M=1
% Three shapes, because the point of a GIVEN grid is that it need not be the one the command would
% have built. An evenly spaced grid is the case where Tauchen's own construction coincides with the
% input, and is the one that can be checked against the no-e_grid call exactly. The other two are
% where a bin-midpoint rule has to actually work: a grid that is not evenly spaced, and one that is
% not centred on the mean.
calibNG.enum=15;
enum=calibNG.enum;
q=3;
calibNG.grid.even=(calibNG.mew+linspace(-q*calibNG.sigma,q*calibNG.sigma,enum))';
% Non-uniform: points clustered near the middle, sparse in the tails. sinh spacing is monotone, so
% the grid is strictly ascending by construction whatever the parameters.
uu=linspace(-1,1,enum)';
calibNG.grid.nonuniform=calibNG.mew+q*calibNG.sigma*sinh(2*uu)/sinh(2);
% Off-centre: the same even grid shifted, so the distribution sits asymmetrically inside it. The
% outermost bins run to +-Inf, so the probabilities must still sum to one.
calibNG.grid.offcentre=calibNG.grid.even+0.8*calibNG.sigma;

%% The multivariate configurations
% M from 1 to 5, because B13 made l_z==5 reachable and nothing has been through it. The sizes are
% chosen so the joint grid stays small: prod(znum) is 15, 49, 125, 256 and 243.
%
% NOT TWO POINTS PER DIMENSION, which is what the first version of this used at M=5. A two-point
% grid puts both points at the edges, +-3sd, so the variance it carries is (3*sd)^2 = 9 times the
% truth whatever the command does - the 2026-09-16 run duly reported a variance error of 13.5
% against a truth of 1.69. That is arithmetic, not a measurement of anything, and it would have
% swamped any real defect at M=5. Three points is the minimum that can place mass at the centre.
calibNG.mv(1).M=1; calibNG.mv(1).znum=15;
calibNG.mv(2).M=2; calibNG.mv(2).znum=[7,7];
calibNG.mv(3).M=3; calibNG.mv(3).znum=[5,5,5];
calibNG.mv(4).M=4; calibNG.mv(4).znum=[4,4,4,4];
calibNG.mv(5).M=5; calibNG.mv(5).znum=[3,3,3,3,3];
for c_c=1:5
    M=calibNG.mv(c_c).M; znum=calibNG.mv(c_c).znum;
    Mew=0.2*(1:M)'-0.3;                       % not zero, and different in each dimension
    sd=0.5+0.2*(0:M-1)';                      % different scales, so a mixed-up dimension shows
    % A correlation matrix that is positive definite for any M: rho^|i-j|, the AR(1) correlation
    % structure, with rho=0.4. Building it this way rather than by hand means the l_z=5 case is as
    % well conditioned as the l_z=2 one.
    R=zeros(M,M);
    for i_c=1:M
        for j_c=1:M
            R(i_c,j_c)=0.4^abs(i_c-j_c);
        end
    end
    Sigma=(sd*sd').*R;
    calibNG.mv(c_c).Mew=Mew; calibNG.mv(c_c).sd=sd; calibNG.mv(c_c).Sigma=Sigma;
    calibNG.mv(c_c).Sigmadiag=diag(diag(Sigma)); % the same variances with the off-diagonals removed
    % The grid, stacked, at +-3 sd of each variable about its own mean
    zg=[];
    for m_c=1:M
        zg=[zg; (Mew(m_c)+linspace(-3*sd(m_c),3*sd(m_c),znum(m_c)))'];
    end
    calibNG.mv(c_c).z_grid=zg;
    fprintf('P9 config M=%i: znum=%s, prod=%i, and the correlation matrix is 0.4^|i-j| \n',M,mat2str(znum),prod(znum))
end

calibNG.tol=10^(-7);
