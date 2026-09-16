function output=DiscP5_crosstests(calib,znums)
% P5: the identities that hold within the VAR(1) block
%
% This is the first block in the bank whose identities depend on a MULTI-DIMENSIONAL state ordering
% being what it is claimed to be, so §1 settles the ordering for both commands before anything else
% relies on it.

fprintf('\n========== P5: cross-tests ========== \n')

output=struct();

%% 1. The state ordering, established directly, and the kron identity
% The two commands index their joint state differently, and every later identity depends on knowing
% which. Rather than infer the convention from the outputs, the two index builders are called on
% distinguishable inputs and asked what they do. That is one line each and it cannot be misread.
fprintf('\n--- which variable moves fastest, in each command''s index builder --- \n')
gv=CreateGridvals([2;3],[10;20;1;2;3],1); % variable 1 has 2 points, variable 2 has 3
fprintf('CreateGridvals (used by discretizeVAR1_Tauchen): first column runs [') 
fprintf('%g ',gv(:,1)); fprintf('] \n')
fprintf('   so variable 1 moves fastest, which is the toolkit convention [T0], this should be one: %i \n',gv(1,1)~=gv(2,1))
ac=allcomb2([10,20,30;1,2,3]);            % two variables, three points each
fprintf('allcomb2 (used by discretizeVAR1_FarmerToda): first column runs [')
fprintf('%g ',ac(:,1)); fprintf('] \n')
fprintf('   so variable 1 moves SLOWEST here, i.e. the LAST variable moves fastest [T0], this should be one: %i \n',ac(1,1)==ac(2,1))
fprintf('   The two commands therefore use OPPOSITE orderings. Only discretizeVAR1_Tauchen matches \n')
fprintf('   CreateGridvals, which is what the rest of the toolkit uses to unpack a joint z grid. \n')
output.T_var1fastest=(gv(1,1)~=gv(2,1));
output.FT_var1fastest=(ac(1,1)~=ac(2,1));

% Now the identity itself, for discretizeVAR1_Tauchen. With Rho and SigmaSq both diagonal the two
% variables are independent, the multivariate normal cdf over a box factorizes into the product of
% its marginal cdf differences, and the grid for each variable is exactly the univariate Tauchen
% grid. So the joint chain must equal the kron of the two univariate Tauchen chains EXACTLY - this
% is a T1 identity, not an accuracy comparison. MATLAB's kron(A,B) runs B fastest, and variable 1 is
% fastest here, so the right form is kron(P2,P1).
fprintf('\n--- diagonal Rho and SigmaSq: the joint chain is a product chain --- \n')
Mew=calib.diag.Mew; Rho=calib.diag.Rho; SigmaSq=calib.diag.SigmaSq;
znumK=9; Tauchen_q=3;
[g1T,P1T]=discretizeAR1_Tauchen(Mew(1),Rho(1,1),sqrt(SigmaSq(1,1)),znumK,Tauchen_q,struct());
[g2T,P2T]=discretizeAR1_Tauchen(Mew(2),Rho(2,2),sqrt(SigmaSq(2,2)),znumK,Tauchen_q,struct());
g1T=gather(g1T); P1T=gather(P1T); g2T=gather(g2T); P2T=gather(P2T);
[zgT,pzT]=discretizeVAR1_Tauchen(Mew,Rho,SigmaSq,znumK,Tauchen_q,struct());
zgT=gather(zgT); pzT=gather(pzT);
d_v1fast=max(abs(pzT-kron(P2T,P1T)),[],'all');
d_v2fast=max(abs(pzT-kron(P1T,P2T)),[],'all');
fprintf('VAR1_Tauchen: its grid blocks equal the univariate Tauchen grids [T1], this should be zero: %2.8e \n',max(abs(zgT-[g1T;g2T])))
fprintf('VAR1_Tauchen: pi_z equals kron(P2,P1), variable 1 fastest [T1], this should be zero: %2.8e \n',d_v1fast)
fprintf('VAR1_Tauchen: and it does NOT equal kron(P1,P2), the other orientation: %2.3e \n',d_v2fast)
fprintf('   (both are printed because a symmetric calibration would make them equal and the check \n')
fprintf('    vacuous; the two variables here differ in both persistence and variance, so they do not) \n')
output.T_kron=[d_v1fast,d_v2fast];

% discretizeVAR1_FarmerToda is NOT compared against the univariate discretizeAR1_FarmerToda here,
% and the reason is worth writing down. Its grid is built in a transformed basis and mapped back
% through C = chol(SigmaSq,'lower')*U with U an orthogonal rotation, so even with a diagonal SigmaSq
% its grid is a rotation of the univariate ones rather than equal to them, and the two chains are
% not comparable state by state. Its product structure is checked where it can be checked exactly -
% the rank-one test on every row of pi_z, in DiscP5_VAR1_FarmerToda - and its ordering is settled
% above by asking allcomb2 directly. What IS worth measuring is how far apart the two commands'
% univariate cases are, which is the next section.

%% 2. M=1: each command against its univariate counterpart
% P2 and P4 both measured the discretizeAR1_FarmerToda vs discretizeVAR1_FarmerToda(M=1) gap; this
% is where it belongs. The grids agree exactly, the transitions do not, because the same absolute
% prior floor kappa=1e-8 is applied to densities built on different scales.
fprintf('\n--- M=1: the VAR commands against their univariate counterparts --- \n')
mew1=calib.full.Mew(1); rho1=calib.full.Rho(1,1); sig1=sqrt(calib.full.SigmaSq(1,1));
for c_c=1:length(znums)
    znum=znums(c_c);
    [gA,pA]=discretizeAR1_Tauchen(mew1,rho1,sig1,znum,Tauchen_q,struct());
    [gB,pB]=discretizeVAR1_Tauchen(mew1,rho1,sig1^2,znum,Tauchen_q,struct());
    gA=gather(gA); pA=gather(pA); gB=gather(gB); pB=gather(pB);
    fprintf('znum=%3i: AR1_Tauchen vs VAR1_Tauchen(M=1) - z_grid [T1] %2.8e, pi_z [T1] %2.8e \n',znum,max(abs(gA(:)-gB(:))),max(abs(pA-pB),[],'all'))
    optsE=struct(); optsE.method='even'; optsE.nSigmas=Tauchen_q; optsE.parallel=1; optsE.verbose=0;
    [gC,pC]=discretizeAR1_FarmerToda(mew1,rho1,sig1,znum,optsE);
    [gD,pD]=discretizeVAR1_FarmerToda(mew1,rho1,sig1^2,znum,optsE);
    gC=gather(gC); pC=gather(pC); gD=gather(gD); pD=gather(pD);
    fprintf('znum=%3i: AR1_FarmerToda vs VAR1_FarmerToda(M=1) - z_grid [T1] %2.8e, pi_z %2.3e (reported; \n',znum,max(abs(gC(:)-gD(:))),max(abs(pC-pD),[],'all'))
    fprintf('   the kappa=1e-8 prior floor is applied on different scales in the two) \n')
end

%% 3. SigmaSq as an M-by-1 column versus diag of it
% Both are documented inputs of discretizeVAR1_Tauchen and must give bit-identical output: the
% command's own first act is to replace the vector by diag of it.
fprintf('\n--- SigmaSq as a column vector versus as a diagonal matrix --- \n')
sv=[0.04;0.09];
[gv,pv]=discretizeVAR1_Tauchen(calib.diag.Mew,calib.diag.Rho,sv,znumK,Tauchen_q,struct());
[gm,pm]=discretizeVAR1_Tauchen(calib.diag.Mew,calib.diag.Rho,diag(sv),znumK,Tauchen_q,struct());
fprintf('z_grid [T0], this should be zero: %2.8e \n',max(abs(gather(gv)-gather(gm))))
fprintf('pi_z [T0], this should be zero: %2.8e \n',max(abs(gather(pv)-gather(pm)),[],'all'))

%% 4. znum and Tauchen_q as scalars versus as vectors (the B26 regression)
% B26, found while scoping this block: the scalar expansion read
%     if isscalar(znum), znum=znum*ones(length(znum),1); end
% and length(scalar) is 1, so it was a no-op and l_z came out as 1 whatever the size of Rho. A
% documented input form silently produced a ONE-variable discretization of a two-variable VAR -
% znum-by-znum instead of znum^2-by-znum^2. Fixed by taking l_z from size(Rho,1). This check is
% cheap and it is the only thing standing between that bug and its return.
fprintf('\n--- scalar versus vector znum and Tauchen_q (the B26 regression) --- \n')
[gs,ps]=discretizeVAR1_Tauchen(calib.full.Mew,calib.full.Rho,calib.full.SigmaSq,znumK,Tauchen_q,struct());
[gvv,pvv]=discretizeVAR1_Tauchen(calib.full.Mew,calib.full.Rho,calib.full.SigmaSq,[znumK;znumK],[Tauchen_q;Tauchen_q],struct());
gs=gather(gs); ps=gather(ps); gvv=gather(gvv); pvv=gather(pvv);
fprintf('scalar znum gives the full znum^M-by-znum^M transition [T0], this should be zero: %i \n',any(size(ps)~=[znumK^2,znumK^2]))
fprintf('scalar znum gives a grid of length sum(znum) [T0], this should be zero: %i \n',any(size(gs)~=[2*znumK,1]))
fprintf('scalar and vector forms agree on z_grid [T0], this should be zero: %2.8e \n',max(abs(gs-gvv)))
fprintf('scalar and vector forms agree on pi_z [T0], this should be zero: %2.8e \n',max(abs(ps-pvv),[],'all'))
% and a genuinely uneven grid, which only the vector form can express
[gu,pu]=discretizeVAR1_Tauchen(calib.full.Mew,calib.full.Rho,calib.full.SigmaSq,[7;11],Tauchen_q,struct());
fprintf('uneven znum=[7;11]: pi_z is 77-by-77 [T0], this should be zero: %i \n',any(size(gather(pu))~=[77,77]))
fprintf('uneven znum=[7;11]: z_grid has length 18 [T0], this should be zero: %i \n',any(size(gather(gu))~=[18,1]))

%% 5. The error paths
fprintf('\n--- error paths --- \n')
try
    discretizeVAR1_Tauchen(calib.full.Mew,[1.01,0;0,0.5],calib.full.SigmaSq,5,Tauchen_q,struct());
    fprintf('non-stationary Rho: ran without error [T0], this should not happen \n')
catch
    fprintf('non-stationary Rho: errors, as it should (this is what B4 fixed - the check used to \n')
    fprintf('   test the eigenvalues without taking absolute values, so a negative unit root passed) \n')
end
try
    discretizeVAR1_Tauchen(calib.full.Mew,calib.full.Rho,[0.04,0.2;0.2,0.09],5,Tauchen_q,struct());
    fprintf('non-positive-definite SigmaSq: ran without error [T0], this should not happen \n')
catch
    fprintf('non-positive-definite SigmaSq: errors, as it should \n')
end
try
    discretizeVAR1_Tauchen(calib.full.Mew,calib.full.Rho,calib.full.SigmaSq,[5;7;9],Tauchen_q,struct());
    fprintf('znum with the wrong number of elements: ran without error [T0], this should not happen \n')
catch
    fprintf('znum with the wrong number of elements: errors, as it should (added with B26) \n')
end

%% 6. The Lyapunov identity, checked on the command's own grid width
% B19 was a wrong unconditional variance. The setup already checks that the truth satisfies its own
% equation; this checks that the command's grid is built from that same quantity, at several
% calibrations rather than just the one the sweep uses, since B19 was wrong by a factor that
% depended on the calibration (5x on one variable, 3x on the other, and 1/sigma at M=1).
fprintf('\n--- the grid width against the Lyapunov solution, at three calibrations --- \n')
cn={'full','diag','three'};
for k_c=1:3
    cc=calib.(cn{k_c});
    if znumK^cc.M>calib.Ncap
        znumL=5;
    else
        znumL=znumK;
    end
    [gL,~]=discretizeVAR1_Tauchen(cc.Mew,cc.Rho,cc.SigmaSq,znumL,Tauchen_q,struct());
    gL=gather(gL);
    hw=zeros(cc.M,1);
    for v_c=1:cc.M
        blk=gL((v_c-1)*znumL+1:v_c*znumL);
        hw(v_c)=(blk(end)-blk(1))/2;
    end
    fprintf('%s (M=%i): |grid half-width/Tauchen_q - sqrt(diag(SigmaSqz))| [T1], this should be zero: %2.8e \n',cn{k_c},cc.M,max(abs(hw/Tauchen_q-cc.sigmaz)))
end

output.done=1;

end
