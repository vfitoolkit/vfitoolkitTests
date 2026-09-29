function output=CoreStationaryGE_InfHorz_PType_GEptype_constraints(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% InfHorz PType: the parameter-constraint sweep, run on a model whose GovBudget condition holds
% CONDITIONAL ON THE PERMANENT TYPE (heteroagentoptions.GEptype), so that Tr is a per-type price.
% Same model as CoreStationaryGE_InfHorz_PType_GEptype: tau*w*N_i=Tr_i, r and tau_c economy-wide.
%
% WHY THIS IS NOT JUST THE EXISTING CONSTRAINT SWEEP AGAIN. A constrained price is transformed to an
% unconstrained one before the optimizer sees it and back again afterwards, and the transform is
% applied by index into the price vector. Under GEptype that vector is no longer one entry per price:
% Tr occupies N_i of them. So this sweep is the one that can catch a transform applied to the wrong
% slots, or applied to only the first of a per-type price's entries. constrainpositive is put on Tr
% precisely because Tr is the per-type price here.
%
% WHICH SOLVERS, AND WHY THESE THREE. fminalgo=1 (fminsearch), 5 (shooting) and 9 (Anderson
% acceleration of the shooting map). Not a full fminalgo sweep: 4 and 8 reach the equilibrium through
% the same objective function as 1 and would only be re-testing that, whereas 5 and 9 never form the
% objective at all - they update each price directly from its paired general eqm condition through
% the howtoupdate rules, so they are a different code path, and they are the path where a per-type
% price has to be expanded into N_i update rules (setupGEnewprice3_shooting does that expansion, and
% errors if a per-type condition is paired with a price that is not itself per-type).
%
% Under 5 and 9 every variant starts from the fminalgo=1 unconstrained answer scaled by 1.05, i.e.
% deliberately off it. That keeps the cost of a shooting solve down without making the test vacuous:
% each variant still has to walk back, through its own constraint transform, to the same place. All
% five variants under a given solver start from the same point, so they must land together.
%
% GE prices: r (capital market), Tr (gov budget, BY PTYPE here), tau_c (consumption-tax budget).
% w is hardcoded from r via the firm FOC, so it is NOT a GE price.

n_p=0;

% Permanent types: two types differing in sigma
N_i=2;
Names_i={'ptype001','ptype002'}; % named, because the per-type price is indexed by name
Params.sigma=[2.2,1.8];       % differs by permanent type
Params.ptypemass=[0.5,0.5];   % mass of each permanent type
PTypeDistParamNames={'ptypemass'};

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi) ...
    ReturnFn_InfHorz(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi);

FnsToEvaluate.K=@(d,aprime,a,z) a;
FnsToEvaluate.N=@(d,aprime,a,z) d*z;
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*z*d+Tr-aprime)/(1+tau_c);

GeneralEqmEqns.CapitalMarket=@(r,K,N,alpha,delta,A) r-(alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta);
GeneralEqmEqns.GovBudget=@(tau,r,N,Tr,alpha,delta,A) tau*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*N-Tr;
GeneralEqmEqns.ConsTax=@(tau_c,C,G) tau_c*C-G;

% Tr as a per-type price, same initial guess for both types
Trstart.ptype001=Params.Tr;
Trstart.ptype002=Params.Tr;
Params.Tr=Trstart;

% The howtoupdate rules for the two update-rule solvers. One row per general eqm condition, exactly
% as without GEptype: the GovBudget row is expanded internally into one rule per permanent type.
%
% THE r FACTOR IS 0.001, NOT THE 0.002 THAT CoreStationaryGE_InfHorz_PType_fminalgo USES, and the
% reason is measured rather than guessed. At 0.002 the 2026-09-23 run showed shooting reach
% CapitalMarket=1.3e-05 within about 200 iterations and then sit in a PERIOD-2 LIMIT CYCLE, the
% condition alternating -1.1e-05, +1.3e-05 every single iteration for the remaining 9,800: constant
% amplitude means the map's multiplier |1-factor*dCapitalMarket/dr| is almost exactly 1, i.e.
% factor*dc/dr is at the stability boundary of 2, so dc/dr is about 1000 for this model (r-MPK is
% steep because K responds strongly to r). The cycle straddles zero with amplitude 2.4e-05, which is
% above toleranceGEcondns=1e-5, so the solve could never exit and ran to maxiter - five variants x
% 10,001 solves was 82%% of that entire 39-hour run. Halving the factor puts the multiplier near zero
% (and it stays well inside 2 even if dc/dr is off by a factor of two either way), so the cycle is
% gone rather than merely smaller. The other two rows are nowhere near the boundary - dc/dTr is -1 and
% dc/dtau_c is C, so their factor*dc/dp is about 0.02 - and their alternation was r's cycle dragging
% them, not an instability of their own.
%
% THE Tr FACTOR IS 0.1, FIVE TIMES THE 0.02 IT WAS, for the reason the FHorz twin raised its own: since
% 2026-09-25 the shooting/Anderson step is taken in the UNCONSTRAINED space, so for a constrained price
% the step in original units is factor*dp/du, and for constrainpositive on Tr with softplus at the
% answer Tr~0.219 that slope is 1-exp(-Tr)=0.196. At 0.02 the three constrainpositive variants were
% therefore stepping five times shorter than the unconstrained one, about 1,800 iterations to converge
% against 350. 0.1 puts the effective step back where it was and stays far inside stability for the
% UNCONSTRAINED variant, where the slope is 1 and dc/dTr is -1, so factor*dc/dTr is 0.1 against a
% boundary of 2. The 'constrain all three' variant is compensated separately in the cc==5 branch below.
howtoupdate={...
    'CapitalMarket','r',0,0.001;   % r_new = r - factor*(r-MPK)
    'GovBudget','Tr',1,0.1;         % Tr_new = Tr + factor*(tau*w*N_i-Tr), one per ptype
    'ConsTax','tau_c',0,0.02};      % tau_c_new = tau_c - factor*(tau_c*C-G)

algonames={'fminalgo=1 (fminsearch)','fminalgo=5 (shooting)  ','fminalgo=9 (Anderson)  '};
varnames={'unconstrained        ','constrainpos default ','constrainpos softplus','constrainpos log     ','constrain all three  '};
p_eqm=cell(3,5);
GEcondns=cell(3,5);

for aa=1:3
    % Solver, and where it starts from
    heteroagentoptionsA=heteroagentoptions;
    heteroagentoptionsA.GEptype={'GovBudget'};
    ParamsA=Params;
    if aa==1
        heteroagentoptionsA.fminalgo=1;
    elseif aa==2
        heteroagentoptionsA.fminalgo=5;
        heteroagentoptionsA.fminalgo5.howtoupdate=howtoupdate;
        heteroagentoptionsA.maxiter=3000; % with the retuned r factor this converges in tens of iterations;
        % 3000 is a backstop so that a solve which does start cycling again costs minutes, not hours
    elseif aa==3
        heteroagentoptionsA.fminalgo=9;
        heteroagentoptionsA.fminalgo9.howtoupdate=howtoupdate;
        heteroagentoptionsA.anderson.maxiter=1e4;
    end
    if aa>1 % start from the fminalgo=1 unconstrained answer, scaled off it
        ParamsA.r=1.05*p_eqm{1,1}.r;
        ParamsA.tau_c=1.05*p_eqm{1,1}.tau_c;
        for ii=1:N_i
            ParamsA.Tr.(Names_i{ii})=1.05*p_eqm{1,1}.Tr.(Names_i{ii});
        end
    end

    for cc=1:5
        heteroagentoptionsC=heteroagentoptionsA;
        if cc==2 % constrainpositive on the per-type price, method left unset (so softplus since 2026-09-12)
            heteroagentoptionsC.constrainpositive={'Tr'};
        elseif cc==3 % the same transform, asked for explicitly: must agree with cc=2 exactly
            heteroagentoptionsC.constrainpositive={'Tr'};
            heteroagentoptionsC.constrainpositivemethod='softplus';
        elseif cc==4 % the other map onto (0,infty): same equilibrium, less accurately on this model
            heteroagentoptionsC.constrainpositive={'Tr'};
            heteroagentoptionsC.constrainpositivemethod='log';
        elseif cc==5 % an economy-wide constrained price on either side of the per-type one
            heteroagentoptionsC.constrain0to1={'r'};
            heteroagentoptionsC.constrainpositive={'Tr'};
            heteroagentoptionsC.constrainAtoB={'tau_c'};
            heteroagentoptionsC.constrainAtoBlimits.tau_c=[0,1];
            % THE r AND tau_c FACTORS ARE COMPENSATED FOR THEIR TRANSFORM SLOPES. The step in original
            % units is factor*dp/du, and here both slopes are small - r(1-r)=0.040 at the answer
            % r=0.0414 under constrain0to1, and (p-a)(b-p)/(b-a)=0.036 at tau_c=0.0378 under
            % constrainAtoB [0,1] - so without this the effective steps are 25x and 28x shorter than
            % the ones the unconstrained variant converges with. On this model that is not merely slow:
            % the FHorz twin, whose r slope is twice as forgiving, still reported maxiter on exactly
            % this variant in the 2026-09-26 run, and r's factor here is already five times smaller.
            % Dividing by the slope puts the effective step back at the base factor, so the variant
            % INHERITS the stability the unconstrained one has - and for r that stability is the whole
            % point of the 0.001 above, whose multiplier |1-factor*dc/dr| sits near zero. The slopes are
            % evaluated at the answer and move about 5% between the 1.05x start and there.
            % Tr needs no bump: the 0.1 in howtoupdate already carries the softplus compensation, and
            % cc=4's log slope is Tr=0.219 against softplus's 0.196, within 12%, so one factor serves both.
            howtoupdateC=howtoupdate;
            howtoupdateC{1,4}=howtoupdate{1,4}/0.040;  % r:     0.001 -> 0.025
            howtoupdateC{3,4}=howtoupdate{3,4}/0.036;  % tau_c: 0.02  -> 0.556
            if aa==2
                heteroagentoptionsC.fminalgo5.howtoupdate=howtoupdateC;
            elseif aa==3
                heteroagentoptionsC.fminalgo9.howtoupdate=howtoupdateC;
            end
        end
        [p_eqm{aa,cc},GEcondns{aa,cc}]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsA, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsC, simoptions, vfoptions);
    end
end

%% Compare
% Each equilibrium is the four numbers [r, tau_c, Tr_ptype001, Tr_ptype002], so a difference is the
% largest gap across all four - the two per-type transfers included, which is the point here.
pvec=zeros(3,5,2+N_i);
for aa=1:3
    for cc=1:5
        pvec(aa,cc,1)=p_eqm{aa,cc}.r;
        pvec(aa,cc,2)=p_eqm{aa,cc}.tau_c;
        for ii=1:N_i
            pvec(aa,cc,2+ii)=p_eqm{aa,cc}.Tr.(Names_i{ii});
        end
    end
end

fprintf('\n=== InfHorz PType GEptype: parameter-constraint invariance ===\n')
for aa=1:3
    fprintf('%s \n',algonames{aa})
    for cc=1:5
        fprintf('  %s r=%.6f tau_c=%.6f Tr=(%.6f, %.6f) \n',varnames{cc},pvec(aa,cc,1),pvec(aa,cc,2),pvec(aa,cc,3),pvec(aa,cc,4))
    end
end
for aa=1:3
    dB=max(abs(pvec(aa,2,:)-pvec(aa,1,:)));
    dB2=max(abs(pvec(aa,3,:)-pvec(aa,1,:)));
    dB3=max(abs(pvec(aa,4,:)-pvec(aa,1,:)));
    dD=max(abs(pvec(aa,5,:)-pvec(aa,1,:)));
    dB2B3=max(abs(pvec(aa,4,:)-pvec(aa,3,:)));
    fprintf('%s constrainpositive on the per-type Tr, default method, this should be near zero: %.3e \n',algonames{aa},dB)
    fprintf('%s constrainpositive on the per-type Tr with softplus, this should be near zero: %.3e \n',algonames{aa},dB2)
    fprintf('%s constrainpositive on the per-type Tr with log, this should be near zero (loose: log is the worse transform on this model): %.3e \n',algonames{aa},dB3)
    fprintf('%s constrain all three, economy-wide and per-type together, this should be near zero: %.3e \n',algonames{aa},dD)
    fprintf('%s log vs softplus, this should NOT be zero: %.3e \n',algonames{aa},dB2B3)
end
% And the two update-rule solvers must reach the same equilibrium as fminsearch did
for aa=2:3
    dU=max(abs(pvec(aa,1,:)-pvec(1,1,:)));
    fprintf('%s unconstrained vs fminalgo=1, this should be near zero: %.3e \n',algonames{aa},dU)
end

output.p_eqm=p_eqm;
output.GEcondns=GEcondns;
output.pvec=pvec;

end
