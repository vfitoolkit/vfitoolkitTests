function output=EZFHorz_bruteforce_oracle(Params,DiscountFactorParamNames)
% Brute-force Epstein-Zin oracle for the main FHorz EZ bank.
%
% Solves a TINY model twice and compares:
%   (a) by explicit nested loops written directly from the Epstein-Zin recursion FORMULAS
%       (gamma, phi, risk aversion appear directly; no ezc coefficients, no vectorisation,
%       no toolkit subcodes — this file shares NOTHING with the toolkit implementation), and
%   (b) through ValueFnIter_Case1_FHorz with the bank's standard EZ vfoptions.
% Every other check in this bank compares two toolkit code paths against each other, so this
% is the only test class that can catch 'toolkit and reference both wrong in the same way'.
%
% Tiny model: n_a=5 (a_grid=linspace(0,2,5)), n_z=2 markov, N_j=3, no d. Ages 1,2 working
% (earnings w*kappa_j*z), age 3 retired (pension); return-fn forms are the bank's three
% EZ ReturnFns (EZReturnFn_cons/positiveUtils/negativeUtils_nod_z_noe_nosemiz),
% reimplemented inline below with plain arithmetic.
%
% Baseline vfoptions convention: vfoptions.EZoneminusbeta defaults to 0, so there is NO
% (1-beta) scaling anywhere in the recursions.
%
% Recursions implemented (written from the formulas; F is the ReturnFn value):
%   Case 1 cons-units (EZutils=0, gamma=ezgamma, phi=ezphi):
%     V = [ F^(1-1/phi) + beta*( E_z[V'^(1-gamma)] )^((1-1/phi)/(1-gamma)) ]^(1/(1-1/phi))
%     terminal: V = [ F^(1-1/phi) ]^(1/(1-1/phi))
%   Case 2 utility-units, positive F (EZutils=1, EZpositiveutility=1, risk=ezrisk):
%     V = F + beta*( E_z[V'^(1-risk)] )^(1/(1-risk)),   terminal: V = F
%   Case 3 utility-units, negative F (EZutils=1, EZpositiveutility=0, risk=ezrisk):
%     V = F - beta*( E_z[(-V')^(1+risk)] )^(1/(1+risk)),   terminal: V = F
% sj + warm-glow variants (survivalprobability + WarmGlowBequestsFn; warm-glow WG(aprime)
% enters the certainty-equivalent with weight (1-sj), and the terminal age uses the
% INSIDE-THE-ROOT warm-glow convention of Kraft-Munk-Weiss 2022):
%   Case 1: V = [ F^(1-1/phi) + beta*( sj*E_z[V'^(1-gamma)] + (1-sj)*WG^(1-gamma) )^((1-1/phi)/(1-gamma)) ]^(1/(1-1/phi))
%     terminal: V = [ F^(1-1/phi) + beta*( (1-sj)*WG^(1-gamma) )^((1-1/phi)/(1-gamma)) ]^(1/(1-1/phi))
%   Case 2: V = F + beta*( sj*E_z[V'^(1-risk)] + (1-sj)*WG^(1-risk) )^(1/(1-risk))
%     terminal: V = F + beta*( (1-sj)*WG^(1-risk) )^(1/(1-risk))
%   Case 3: V = F - beta*( sj*E_z[(-V')^(1+risk)] + (1-sj)*(-WG)^(1+risk) )^(1/(1+risk))
%     terminal: V = F - beta*( (1-sj)*(-WG)^(1+risk) )^(1/(1+risk))
% The warm-glow forms are the bank's EZWarmGlowFn_cons/positiveUtils/negativeUtils,
% reimplemented inline (cons: strictly positive; posU: strictly positive; negU: strictly
% negative — so the sign conventions above are well-defined).
%
% Comparisons: V by max abs diff (expect ~1e-14 roundoff: the brute force orders the
% floating-point operations differently), Policy exactly (both sides break exact ties by
% the first index, and F does not tie across aprime in this model).
%
% All Params changes are made on a local copy P (the caller's Params is untouched).

%% Tiny model
n_a=5;
a_grid=linspace(0,2,n_a)'; % column
n_z=2;
z_grid=[0.8;1.2];
pi_z=[0.9,0.1;0.2,0.8];
N_j=3;

P=Params; % local copy: age-dependent params must be length N_j=3
P.agej=1:1:N_j;
P.Jr=3;             % ages 1,2 working; age 3 retired (pension)
P.kappa_j=[0.5,1,0];
P.sj=[0.9,0.7,0];   % declining survival, sj(N_j)=0 (only used in the sj+warm-glow variants)
% beta, r, w, pension, ezgamma, ezphi, ezrisk, ezsigma, wg1, wg2, wg3 are taken from Params

% Bank ReturnFns (toolkit side; the brute force below reimplements these inline)
ReturnFn_cons=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) EZReturnFn_cons_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,agej,Jr,pension);
ReturnFn_posU=@(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension) EZReturnFn_positiveUtils_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_negU=@(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension) EZReturnFn_negativeUtils_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension);
% Bank warm-glow fns (toolkit side)
WGFn_cons=@(aprime,wg1,wg2) EZWarmGlowFn_cons(aprime,wg1,wg2);
WGFn_posU=@(aprime,wg1,wg2,wg3) EZWarmGlowFn_positiveUtils(aprime,wg1,wg2,wg3);
WGFn_negU=@(aprime,wg1,wg2,wg3) EZWarmGlowFn_negativeUtils(aprime,wg1,wg2,wg3);

% Shorthands for the brute force (plain CPU doubles throughout)
beta=P.beta;
r=P.r;
w=P.w;
pension=P.pension;
kappa_j=P.kappa_j;
Jr=P.Jr;
agej=P.agej;
sj=P.sj;
ezgamma=P.ezgamma;
ezphi=P.ezphi;
ezrisk=P.ezrisk;
ezsigma=P.ezsigma;
wg1=P.wg1;
wg2=P.wg2;
wg3=P.wg3;

fprintf('Brute-force EZ oracle (tiny model: n_a=%i, n_z=%i, N_j=%i, no d) \n',n_a,n_z,N_j)

%% Case 1: consumption-units — brute force from the formula
rho=1-1/ezphi;    % exponent on the composite consumption good
alpha=1-ezgamma;  % exponent inside the certainty-equivalent
Vb=zeros(n_a,n_z,N_j);
Pb=ones(n_a,n_z,N_j);
for jj=N_j:-1:1
    for i_z=1:n_z
        for i_a=1:n_a
            bestV=-Inf;
            besti=1;
            for i_ap=1:n_a
                if agej(jj)<Jr
                    c=(1+r)*a_grid(i_a)+w*kappa_j(jj)*z_grid(i_z)-a_grid(i_ap);
                else
                    c=(1+r)*a_grid(i_a)+pension-a_grid(i_ap);
                end
                if c>0
                    if jj==N_j
                        Vcand=(c^rho)^(1/rho);
                    else
                        EVpow=0;
                        for i_zp=1:n_z
                            EVpow=EVpow+pi_z(i_z,i_zp)*Vb(i_ap,i_zp,jj+1)^alpha;
                        end
                        Vcand=(c^rho+beta*EVpow^(rho/alpha))^(1/rho);
                    end
                    if Vcand>bestV
                        bestV=Vcand;
                        besti=i_ap;
                    end
                end
            end
            Vb(i_a,i_z,jj)=bestV;
            Pb(i_a,i_z,jj)=besti;
        end
    end
end

% Toolkit solve with the bank's standard cons-units EZ vfoptions
vfoptions1=struct();
vfoptions1.exoticpreferences='EpsteinZin';
vfoptions1.EZutils=0;
vfoptions1.EZriskaversion='ezgamma';
vfoptions1.EZeis='ezphi';
[V1,Policy1]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn_cons,P,DiscountFactorParamNames,[],vfoptions1);
V1=reshape(gather(V1),[n_a,n_z,N_j]);
Policy1=reshape(gather(Policy1),[n_a,n_z,N_j]);
fprintf('Brute-force oracle vs toolkit, V [EZ cons-units], should be zero up to ~1e-14 roundoff: %g \n',max(abs(V1(:)-Vb(:))))
fprintf('Brute-force oracle vs toolkit, Policy [EZ cons-units], this should be zero: %2.8f \n',max(abs(Policy1(:)-Pb(:))))

%% Case 2: utility-units, POSITIVE utility fn — brute force from the formula
riskexp=1-ezrisk;
Vb=zeros(n_a,n_z,N_j);
Pb=ones(n_a,n_z,N_j);
for jj=N_j:-1:1
    for i_z=1:n_z
        for i_a=1:n_a
            bestV=-Inf;
            besti=1;
            for i_ap=1:n_a
                if agej(jj)<Jr
                    c=(1+r)*a_grid(i_a)+w*kappa_j(jj)*z_grid(i_z)-a_grid(i_ap);
                else
                    c=(1+r)*a_grid(i_a)+pension-a_grid(i_ap);
                end
                if c>0
                    F=((1+c)^(1-ezsigma)-1)/(1-ezsigma); % strictly positive
                    if jj==N_j
                        Vcand=F;
                    else
                        EVpow=0;
                        for i_zp=1:n_z
                            EVpow=EVpow+pi_z(i_z,i_zp)*Vb(i_ap,i_zp,jj+1)^riskexp;
                        end
                        Vcand=F+beta*EVpow^(1/riskexp);
                    end
                    if Vcand>bestV
                        bestV=Vcand;
                        besti=i_ap;
                    end
                end
            end
            Vb(i_a,i_z,jj)=bestV;
            Pb(i_a,i_z,jj)=besti;
        end
    end
end

vfoptions2=struct();
vfoptions2.exoticpreferences='EpsteinZin';
vfoptions2.EZutils=1;
vfoptions2.EZpositiveutility=1;
vfoptions2.EZriskaversion='ezrisk';
[V2,Policy2]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn_posU,P,DiscountFactorParamNames,[],vfoptions2);
V2=reshape(gather(V2),[n_a,n_z,N_j]);
Policy2=reshape(gather(Policy2),[n_a,n_z,N_j]);
fprintf('Brute-force oracle vs toolkit, V [EZ positive utils], should be zero up to ~1e-14 roundoff: %g \n',max(abs(V2(:)-Vb(:))))
fprintf('Brute-force oracle vs toolkit, Policy [EZ positive utils], this should be zero: %2.8f \n',max(abs(Policy2(:)-Pb(:))))

%% Case 3: utility-units, NEGATIVE utility fn — brute force from the formula
q=1+ezrisk;
Vb=zeros(n_a,n_z,N_j);
Pb=ones(n_a,n_z,N_j);
for jj=N_j:-1:1
    for i_z=1:n_z
        for i_a=1:n_a
            bestV=-Inf;
            besti=1;
            for i_ap=1:n_a
                if agej(jj)<Jr
                    c=(1+r)*a_grid(i_a)+w*kappa_j(jj)*z_grid(i_z)-a_grid(i_ap);
                else
                    c=(1+r)*a_grid(i_a)+pension-a_grid(i_ap);
                end
                if c>0
                    F=(c^(1-ezsigma))/(1-ezsigma); % strictly negative
                    if jj==N_j
                        Vcand=F;
                    else
                        EVpow=0;
                        for i_zp=1:n_z
                            EVpow=EVpow+pi_z(i_z,i_zp)*(-Vb(i_ap,i_zp,jj+1))^q;
                        end
                        Vcand=F-beta*EVpow^(1/q);
                    end
                    if Vcand>bestV
                        bestV=Vcand;
                        besti=i_ap;
                    end
                end
            end
            Vb(i_a,i_z,jj)=bestV;
            Pb(i_a,i_z,jj)=besti;
        end
    end
end

vfoptions3=struct();
vfoptions3.exoticpreferences='EpsteinZin';
vfoptions3.EZutils=1;
vfoptions3.EZpositiveutility=0;
vfoptions3.EZriskaversion='ezrisk';
[V3,Policy3]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn_negU,P,DiscountFactorParamNames,[],vfoptions3);
V3=reshape(gather(V3),[n_a,n_z,N_j]);
Policy3=reshape(gather(Policy3),[n_a,n_z,N_j]);
fprintf('Brute-force oracle vs toolkit, V [EZ negative utils], should be zero up to ~1e-14 roundoff: %g \n',max(abs(V3(:)-Vb(:))))
fprintf('Brute-force oracle vs toolkit, Policy [EZ negative utils], this should be zero: %2.8f \n',max(abs(Policy3(:)-Pb(:))))

%% sj + warm-glow, Case 1: consumption-units — brute force from the formula
% Independently validates the inside-the-root terminal warm-glow convention (Kraft-Munk-Weiss)
rho=1-1/ezphi;
alpha=1-ezgamma;
Vb=zeros(n_a,n_z,N_j);
Pb=ones(n_a,n_z,N_j);
for jj=N_j:-1:1
    for i_z=1:n_z
        for i_a=1:n_a
            bestV=-Inf;
            besti=1;
            for i_ap=1:n_a
                if agej(jj)<Jr
                    c=(1+r)*a_grid(i_a)+w*kappa_j(jj)*z_grid(i_z)-a_grid(i_ap);
                else
                    c=(1+r)*a_grid(i_a)+pension-a_grid(i_ap);
                end
                if c>0
                    WG=wg1*(1+a_grid(i_ap)/wg2); % cons-units warm-glow, strictly positive
                    if jj==N_j
                        Vcand=(c^rho+beta*((1-sj(N_j))*WG^alpha)^(rho/alpha))^(1/rho);
                    else
                        EVpow=0;
                        for i_zp=1:n_z
                            EVpow=EVpow+pi_z(i_z,i_zp)*Vb(i_ap,i_zp,jj+1)^alpha;
                        end
                        Vcand=(c^rho+beta*(sj(jj)*EVpow+(1-sj(jj))*WG^alpha)^(rho/alpha))^(1/rho);
                    end
                    if Vcand>bestV
                        bestV=Vcand;
                        besti=i_ap;
                    end
                end
            end
            Vb(i_a,i_z,jj)=bestV;
            Pb(i_a,i_z,jj)=besti;
        end
    end
end

vfoptions1wg=vfoptions1;
vfoptions1wg.survivalprobability='sj';
vfoptions1wg.WarmGlowBequestsFn=WGFn_cons;
[V4,Policy4]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn_cons,P,DiscountFactorParamNames,[],vfoptions1wg);
V4=reshape(gather(V4),[n_a,n_z,N_j]);
Policy4=reshape(gather(Policy4),[n_a,n_z,N_j]);
fprintf('Brute-force oracle vs toolkit, sj+warm-glow V [EZ cons-units], should be zero up to ~1e-14 roundoff: %g \n',max(abs(V4(:)-Vb(:))))
fprintf('Brute-force oracle vs toolkit, sj+warm-glow Policy [EZ cons-units], this should be zero: %2.8f \n',max(abs(Policy4(:)-Pb(:))))

%% sj + warm-glow, Case 2: utility-units, POSITIVE utility fn — brute force from the formula
riskexp=1-ezrisk;
Vb=zeros(n_a,n_z,N_j);
Pb=ones(n_a,n_z,N_j);
for jj=N_j:-1:1
    for i_z=1:n_z
        for i_a=1:n_a
            bestV=-Inf;
            besti=1;
            for i_ap=1:n_a
                if agej(jj)<Jr
                    c=(1+r)*a_grid(i_a)+w*kappa_j(jj)*z_grid(i_z)-a_grid(i_ap);
                else
                    c=(1+r)*a_grid(i_a)+pension-a_grid(i_ap);
                end
                if c>0
                    F=((1+c)^(1-ezsigma)-1)/(1-ezsigma); % strictly positive
                    WG=wg1*((2+a_grid(i_ap)/wg2)^(1-wg3)-1)/(1-wg3); % strictly positive
                    if jj==N_j
                        Vcand=F+beta*((1-sj(N_j))*WG^riskexp)^(1/riskexp);
                    else
                        EVpow=0;
                        for i_zp=1:n_z
                            EVpow=EVpow+pi_z(i_z,i_zp)*Vb(i_ap,i_zp,jj+1)^riskexp;
                        end
                        Vcand=F+beta*(sj(jj)*EVpow+(1-sj(jj))*WG^riskexp)^(1/riskexp);
                    end
                    if Vcand>bestV
                        bestV=Vcand;
                        besti=i_ap;
                    end
                end
            end
            Vb(i_a,i_z,jj)=bestV;
            Pb(i_a,i_z,jj)=besti;
        end
    end
end

vfoptions2wg=vfoptions2;
vfoptions2wg.survivalprobability='sj';
vfoptions2wg.WarmGlowBequestsFn=WGFn_posU;
[V5,Policy5]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn_posU,P,DiscountFactorParamNames,[],vfoptions2wg);
V5=reshape(gather(V5),[n_a,n_z,N_j]);
Policy5=reshape(gather(Policy5),[n_a,n_z,N_j]);
fprintf('Brute-force oracle vs toolkit, sj+warm-glow V [EZ positive utils], should be zero up to ~1e-14 roundoff: %g \n',max(abs(V5(:)-Vb(:))))
fprintf('Brute-force oracle vs toolkit, sj+warm-glow Policy [EZ positive utils], this should be zero: %2.8f \n',max(abs(Policy5(:)-Pb(:))))

%% sj + warm-glow, Case 3: utility-units, NEGATIVE utility fn — brute force from the formula
q=1+ezrisk;
Vb=zeros(n_a,n_z,N_j);
Pb=ones(n_a,n_z,N_j);
for jj=N_j:-1:1
    for i_z=1:n_z
        for i_a=1:n_a
            bestV=-Inf;
            besti=1;
            for i_ap=1:n_a
                if agej(jj)<Jr
                    c=(1+r)*a_grid(i_a)+w*kappa_j(jj)*z_grid(i_z)-a_grid(i_ap);
                else
                    c=(1+r)*a_grid(i_a)+pension-a_grid(i_ap);
                end
                if c>0
                    F=(c^(1-ezsigma))/(1-ezsigma); % strictly negative
                    WG=wg1*((1+a_grid(i_ap)/wg2)^(1-wg3))/(1-wg3); % strictly negative
                    if jj==N_j
                        Vcand=F-beta*((1-sj(N_j))*(-WG)^q)^(1/q);
                    else
                        EVpow=0;
                        for i_zp=1:n_z
                            EVpow=EVpow+pi_z(i_z,i_zp)*(-Vb(i_ap,i_zp,jj+1))^q;
                        end
                        Vcand=F-beta*(sj(jj)*EVpow+(1-sj(jj))*(-WG)^q)^(1/q);
                    end
                    if Vcand>bestV
                        bestV=Vcand;
                        besti=i_ap;
                    end
                end
            end
            Vb(i_a,i_z,jj)=bestV;
            Pb(i_a,i_z,jj)=besti;
        end
    end
end

vfoptions3wg=vfoptions3;
vfoptions3wg.survivalprobability='sj';
vfoptions3wg.WarmGlowBequestsFn=WGFn_negU;
[V6,Policy6]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn_negU,P,DiscountFactorParamNames,[],vfoptions3wg);
V6=reshape(gather(V6),[n_a,n_z,N_j]);
Policy6=reshape(gather(Policy6),[n_a,n_z,N_j]);
fprintf('Brute-force oracle vs toolkit, sj+warm-glow V [EZ negative utils], should be zero up to ~1e-14 roundoff: %g \n',max(abs(V6(:)-Vb(:))))
fprintf('Brute-force oracle vs toolkit, sj+warm-glow Policy [EZ negative utils], this should be zero: %2.8f \n',max(abs(Policy6(:)-Pb(:))))

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
