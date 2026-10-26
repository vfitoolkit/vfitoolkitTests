function output=EZInfHorz_bruteforce_oracle(Params,DiscountFactorParamNames)
% Brute-force Epstein-Zin oracle for the InfHorz EZ bank (the InfHorz twin of EZFHorz_bruteforce_oracle).
%
% Solves a TINY model twice and compares:
%   (a) by explicit nested loops written directly from the Epstein-Zin recursion FORMULAS
%       (gamma, phi, risk appear directly; no ezc coefficients, no vectorisation, no toolkit
%       subcodes: this file shares NOTHING with the toolkit implementation), iterated to a fixed
%       point from the 'consume everything' value, and
%   (b) through ValueFnIter_InfHorz with the bank's standard EZ vfoptions (which starts from V=0
%       and uses Howards; a different path to the same fixed point).
% Every other check in this bank compares two toolkit code paths against each other, or the
% toolkit against an identity evaluated at the toolkit's own answer, so this is the only test
% class that can catch 'toolkit and reference both wrong in the same way'.
%
% Tiny model: n_a=5 (a_grid=linspace(0,2,5)), n_z=2 markov, no d, c=(1+r)a+w*z-aprime; the
% return-fn forms are the bank's three EZ ReturnFns (EZReturnFn_cons/positiveUtils/negativeUtils_nod_z),
% reimplemented inline below with plain arithmetic. r, w, beta and the EZ parameters are
% taken from Params. vfoptions.EZoneminusbeta defaults to 0, so there is NO (1-beta) scaling.
%
% Recursions implemented (F is the ReturnFn value):
%   Case 1 cons-units (EZutils=0, gamma=ezgamma, phi=ezphi):
%     V = [ F^(1-1/phi) + beta*( E_z[V'^(1-gamma)] )^((1-1/phi)/(1-gamma)) ]^(1/(1-1/phi))
%   Case 2 utility-units, positive F (EZutils=1, EZpositiveutility=1, risk=ezrisk):
%     V = F + beta*( E_z[V'^(1-risk)] )^(1/(1-risk))
%   Case 3 utility-units, negative F (EZutils=1, EZpositiveutility=0, risk=ezrisk):
%     V = F - beta*( E_z[(-V')^(1+risk)] )^(1/(1+risk))
%
% Comparisons: V by max abs diff (both sides are fixed points to tight tolerances: the brute
% force to 1e-13 on the V change, the toolkit to vfoptions.tolerance=1e-12); Policy exactly
% (no ties in this model), AND the brute-force objective evaluated at the toolkit's policy
% against the brute-force V, which is the tie-proof form of the same claim: if the toolkit
% policy is optimal for the brute-force recursion, the brute-force objective at it equals the
% brute-force V to convergence tolerance, whatever index a tie resolved to.

%% Tiny model
n_a=5;
a_grid=linspace(0,2,n_a)'; % column
n_z=2;
z_grid=[0.8;1.2];
pi_z=[0.9,0.1;0.2,0.8];

P=Params; % local copy (nothing is changed, but keeps the caller's Params out of the loops)

% Bank ReturnFns (toolkit side; the brute force below reimplements these inline)
ReturnFn_cons=@(aprime,a,z,r,w) EZReturnFn_cons_nod_z(aprime,a,z,r,w);
ReturnFn_posU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_positiveUtils_nod_z(aprime,a,z,r,w,ezsigma);
ReturnFn_negU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_negativeUtils_nod_z(aprime,a,z,r,w,ezsigma);

% Shorthands for the brute force (plain CPU doubles throughout)
beta=P.beta;
r=P.r;
w=P.w;
ezgamma=P.ezgamma;
ezphi=P.ezphi;
ezrisk=P.ezrisk;
ezsigma=P.ezsigma;

brutetol=10^(-13);
brutemaxiter=10^5;
vftol=10^(-12); % the toolkit solves, tighter than the rest of the bank so that the two fixed points can be compared at the identity bar

fprintf('Brute-force EZ oracle (tiny InfHorz model: n_a=%i, n_z=%i, no d) \n',n_a,n_z)

%% Case 1: consumption-units, brute force from the formula
rho=1-1/ezphi;    % exponent on the composite consumption good
alpha=1-ezgamma;  % exponent inside the certainty-equivalent
% Start from the 'consume everything' value (aprime=0, V=F=c): a feasible, positive starting point
Vb=zeros(n_a,n_z);
for i_z=1:n_z
    for i_a=1:n_a
        Vb(i_a,i_z)=(1+r)*a_grid(i_a)+w*z_grid(i_z)-a_grid(1);
    end
end
Pb=ones(n_a,n_z);
bdist=Inf;
biter=0;
while bdist>brutetol && biter<brutemaxiter
    Vbold=Vb;
    for i_z=1:n_z
        for i_a=1:n_a
            bestV=-Inf;
            besti=1;
            for i_ap=1:n_a
                c=(1+r)*a_grid(i_a)+w*z_grid(i_z)-a_grid(i_ap);
                if c>0
                    EVpow=0;
                    for i_zp=1:n_z
                        EVpow=EVpow+pi_z(i_z,i_zp)*Vbold(i_ap,i_zp)^alpha;
                    end
                    Vcand=(c^rho+beta*EVpow^(rho/alpha))^(1/rho);
                    if Vcand>bestV
                        bestV=Vcand;
                        besti=i_ap;
                    end
                end
            end
            Vb(i_a,i_z)=bestV;
            Pb(i_a,i_z)=besti;
        end
    end
    bdist=max(abs(Vb(:)-Vbold(:)));
    biter=biter+1;
end
fprintf('Brute force [EZ cons-units] converged in %i iterations (last change %g); converged, this should be one: %i \n',biter,bdist,bdist<=brutetol)
% Toolkit solve with the bank's standard cons-units EZ vfoptions
vfoptions1=struct();
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=0;
vfoptions1.EZriskaversion='ezgamma';
vfoptions1.EZeis='ezphi';
vfoptions1.tolerance=vftol;
[V1,Policy1]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,P,DiscountFactorParamNames,[],vfoptions1);
V1=reshape(gather(V1),[n_a,n_z]);
Policy1=reshape(gather(Policy1),[n_a,n_z]);
fprintf('Brute-force oracle vs toolkit, V [EZ cons-units], this should be zero to convergence tolerance: %2.8e \n',max(abs(V1(:)-Vb(:))))
fprintf('Brute-force oracle vs toolkit, Policy [EZ cons-units], this should be zero: %.3e \n',max(abs(Policy1(:)-Pb(:))))
% The brute-force objective at the toolkit's policy (tie-proof form of the Policy check)
RHSatP=zeros(n_a,n_z);
for i_z=1:n_z
    for i_a=1:n_a
        i_ap=Policy1(i_a,i_z);
        c=(1+r)*a_grid(i_a)+w*z_grid(i_z)-a_grid(i_ap);
        EVpow=0;
        for i_zp=1:n_z
            EVpow=EVpow+pi_z(i_z,i_zp)*Vb(i_ap,i_zp)^alpha;
        end
        RHSatP(i_a,i_z)=(c^rho+beta*EVpow^(rho/alpha))^(1/rho);
    end
end
fprintf('Brute-force objective at the toolkit policy vs brute-force V [EZ cons-units], this should be zero to convergence tolerance: %2.8e \n',max(abs(RHSatP(:)-Vb(:))))

%% Case 2: utility-units, POSITIVE utility fn, brute force from the formula
riskexp=1-ezrisk;
Vb=zeros(n_a,n_z);
for i_z=1:n_z
    for i_a=1:n_a
        c=(1+r)*a_grid(i_a)+w*z_grid(i_z)-a_grid(1);
        Vb(i_a,i_z)=((1+c)^(1-ezsigma)-1)/(1-ezsigma); % strictly positive
    end
end
Pb=ones(n_a,n_z);
bdist=Inf;
biter=0;
while bdist>brutetol && biter<brutemaxiter
    Vbold=Vb;
    for i_z=1:n_z
        for i_a=1:n_a
            bestV=-Inf;
            besti=1;
            for i_ap=1:n_a
                c=(1+r)*a_grid(i_a)+w*z_grid(i_z)-a_grid(i_ap);
                if c>0
                    F=((1+c)^(1-ezsigma)-1)/(1-ezsigma); % strictly positive
                    EVpow=0;
                    for i_zp=1:n_z
                        EVpow=EVpow+pi_z(i_z,i_zp)*Vbold(i_ap,i_zp)^riskexp;
                    end
                    Vcand=F+beta*EVpow^(1/riskexp);
                    if Vcand>bestV
                        bestV=Vcand;
                        besti=i_ap;
                    end
                end
            end
            Vb(i_a,i_z)=bestV;
            Pb(i_a,i_z)=besti;
        end
    end
    bdist=max(abs(Vb(:)-Vbold(:)));
    biter=biter+1;
end
fprintf('Brute force [EZ positive utils] converged in %i iterations (last change %g); converged, this should be one: %i \n',biter,bdist,bdist<=brutetol)
vfoptions2=struct();
vfoptions2.exoticpreferences='EpsteinZin';
vfoptions2.EZutils=1;
vfoptions2.EZpositiveutility=1;
vfoptions2.EZriskaversion='ezrisk';
vfoptions2.tolerance=vftol;
[V2,Policy2]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_posU,P,DiscountFactorParamNames,[],vfoptions2);
V2=reshape(gather(V2),[n_a,n_z]);
Policy2=reshape(gather(Policy2),[n_a,n_z]);
fprintf('Brute-force oracle vs toolkit, V [EZ positive utils], this should be zero to convergence tolerance: %2.8e \n',max(abs(V2(:)-Vb(:))))
fprintf('Brute-force oracle vs toolkit, Policy [EZ positive utils], this should be zero: %.3e \n',max(abs(Policy2(:)-Pb(:))))
RHSatP=zeros(n_a,n_z);
for i_z=1:n_z
    for i_a=1:n_a
        i_ap=Policy2(i_a,i_z);
        c=(1+r)*a_grid(i_a)+w*z_grid(i_z)-a_grid(i_ap);
        F=((1+c)^(1-ezsigma)-1)/(1-ezsigma);
        EVpow=0;
        for i_zp=1:n_z
            EVpow=EVpow+pi_z(i_z,i_zp)*Vb(i_ap,i_zp)^riskexp;
        end
        RHSatP(i_a,i_z)=F+beta*EVpow^(1/riskexp);
    end
end
fprintf('Brute-force objective at the toolkit policy vs brute-force V [EZ positive utils], this should be zero to convergence tolerance: %2.8e \n',max(abs(RHSatP(:)-Vb(:))))

%% Case 3: utility-units, NEGATIVE utility fn, brute force from the formula
riskexp=1+ezrisk;
Vb=zeros(n_a,n_z);
for i_z=1:n_z
    for i_a=1:n_a
        c=(1+r)*a_grid(i_a)+w*z_grid(i_z)-a_grid(1);
        Vb(i_a,i_z)=(c^(1-ezsigma))/(1-ezsigma); % strictly negative
    end
end
Pb=ones(n_a,n_z);
bdist=Inf;
biter=0;
while bdist>brutetol && biter<brutemaxiter
    Vbold=Vb;
    for i_z=1:n_z
        for i_a=1:n_a
            bestV=-Inf;
            besti=1;
            for i_ap=1:n_a
                c=(1+r)*a_grid(i_a)+w*z_grid(i_z)-a_grid(i_ap);
                if c>0
                    F=(c^(1-ezsigma))/(1-ezsigma); % strictly negative
                    EVpow=0;
                    for i_zp=1:n_z
                        EVpow=EVpow+pi_z(i_z,i_zp)*(-Vbold(i_ap,i_zp))^riskexp;
                    end
                    Vcand=F-beta*EVpow^(1/riskexp);
                    if Vcand>bestV
                        bestV=Vcand;
                        besti=i_ap;
                    end
                end
            end
            Vb(i_a,i_z)=bestV;
            Pb(i_a,i_z)=besti;
        end
    end
    bdist=max(abs(Vb(:)-Vbold(:)));
    biter=biter+1;
end
fprintf('Brute force [EZ negative utils] converged in %i iterations (last change %g); converged, this should be one: %i \n',biter,bdist,bdist<=brutetol)
vfoptions3=struct();
vfoptions3.exoticpreferences='EpsteinZin';
vfoptions3.EZutils=1;
vfoptions3.EZpositiveutility=0;
vfoptions3.EZriskaversion='ezrisk';
vfoptions3.tolerance=vftol;
[V3,Policy3]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_negU,P,DiscountFactorParamNames,[],vfoptions3);
V3=reshape(gather(V3),[n_a,n_z]);
Policy3=reshape(gather(Policy3),[n_a,n_z]);
fprintf('Brute-force oracle vs toolkit, V [EZ negative utils], this should be zero to convergence tolerance: %2.8e \n',max(abs(V3(:)-Vb(:))))
fprintf('Brute-force oracle vs toolkit, Policy [EZ negative utils], this should be zero: %.3e \n',max(abs(Policy3(:)-Pb(:))))
RHSatP=zeros(n_a,n_z);
for i_z=1:n_z
    for i_a=1:n_a
        i_ap=Policy3(i_a,i_z);
        c=(1+r)*a_grid(i_a)+w*z_grid(i_z)-a_grid(i_ap);
        F=(c^(1-ezsigma))/(1-ezsigma);
        EVpow=0;
        for i_zp=1:n_z
            EVpow=EVpow+pi_z(i_z,i_zp)*(-Vb(i_ap,i_zp))^riskexp;
        end
        RHSatP(i_a,i_z)=F-beta*EVpow^(1/riskexp);
    end
end
fprintf('Brute-force objective at the toolkit policy vs brute-force V [EZ negative utils], this should be zero to convergence tolerance: %2.8e \n',max(abs(RHSatP(:)-Vb(:))))

%% Anti-vacuity: the three cases must not share a policy with each other by construction
% (they are different preferences over the same budget set; on five grid points they may still
% agree, so this is a diagnostic rather than a check)
fprintf('DIAGNOSTIC: states where the cons-units and positive-utils policies differ: %i of %i; cons-units vs negative-utils: %i of %i \n',sum(Policy1(:)~=Policy2(:)),n_a*n_z,sum(Policy1(:)~=Policy3(:)),n_a*n_z)

%%
output=struct();

end
