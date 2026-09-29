function output=CoreStationaryGE_FHorz_PType_GEptype_constraints(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% FHorz PType: the parameter-constraint sweep, run on a model whose GovBudget condition holds
% CONDITIONAL ON THE PERMANENT TYPE (heteroagentoptions.GEptype), so that Tr is a per-type price.
% Same model as CoreStationaryGE_FHorz_PType_GEptype: tau*w*N_i=Tr_i, with r and tau_c economy-wide.
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
vfoptions.divideandconquer=1; % finite-horizon tests use divide-and-conquer

% Permanent types: two types differing in sigma
N_i=2;
Names_i={'ptype001','ptype002'}; % named, because the per-type price is indexed by name
Params.sigma=[2.2,1.8];       % differs by permanent type
Params.ptypemass=[0.5,0.5];   % mass of each permanent type
PTypeDistParamNames={'ptypemass'};

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi) ...
    ReturnFn_FHorz(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi);

FnsToEvaluate.K=@(d,aprime,a,z) a;
FnsToEvaluate.N=@(d,aprime,a,z,kappa_j) kappa_j*z*d;
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*kappa_j*z*d+Tr-aprime)/(1+tau_c);

GeneralEqmEqns.CapitalMarket=@(r,K,N,alpha,delta,A) r-(alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta);
GeneralEqmEqns.GovBudget=@(tau,r,N,Tr,alpha,delta,A) tau*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*N-Tr;
GeneralEqmEqns.ConsTax=@(tau_c,C,G) tau_c*C-G;

% Tr as a per-type price, same initial guess for both types
Trstart.ptype001=Params.Tr;
Trstart.ptype002=Params.Tr;
Params.Tr=Trstart;

% The howtoupdate rules for the two update-rule solvers. One row per general eqm condition, exactly
% as without GEptype: the GovBudget row is expanded internally into one rule per permanent type.
% THE Tr FACTOR IS 0.3, SIX TIMES THE 0.05 THE OTHER SUBCODES USE, because since 2026-09-25 the
% shooting/Anderson step is taken in the UNCONSTRAINED space (which is what makes the constraints in
% this very sweep do anything). For a constrained price the step in original units is then
% factor*dp/du, and for constrainpositive on Tr with softplus at Tr~0.184 the slope dp/du is
% sigmoid(u)~0.164 - so the old factor would have taken about six times as many iterations on the
% three constrainpositive variants. 0.3 puts the effective step back where it was, and is still far
% inside stability for the UNCONSTRAINED variant, where the slope is 1 and dc/dTr is -1: factor*dc/dTr
% is 0.3 against a boundary of 2.
%
% The 'constrain all three' variant has r and tau_c constrained too, with smaller slopes still, and it
% is compensated SEPARATELY down in the cc==5 branch rather than by this shared factor - one factor
% cannot serve both that variant and the unconstrained one. Until 2026-09-27 it was not compensated at
% all, and the 2026-09-26 run showed the cost: it was the only variant to report maxiter.
howtoupdate={...
    'CapitalMarket','r',0,0.005;   % r_new = r - factor*(r-MPK)
    'GovBudget','Tr',1,0.3;         % Tr_new = Tr + factor*(tau*w*N_i-Tr), one per ptype
    'ConsTax','tau_c',0,0.05};      % tau_c_new = tau_c - factor*(tau_c*C-G)

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
        heteroagentoptionsA.maxiter=3000; % backstop, as in the InfHorz twin: a solve that cycles costs
        % minutes rather than the 11 hours that 1e4 would at about 4 seconds a solve here
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
            % THE r AND tau_c FACTORS ARE COMPENSATED FOR THEIR TRANSFORM SLOPES, which the comment on
            % howtoupdate above used to say was deliberately not done. The 2026-09-26 run showed what
            % that costs: this was the only variant in the part to report maxiter under fminalgo=5
            % (3000 iterations, not converged), and under fminalgo=9 it took 4,100+ iterations against
            % about 230 for every other variant. The step in original units is factor*dp/du and both
            % slopes are small - r(1-r)=0.077 at the answer r=0.084 under constrain0to1, and
            % (p-a)(b-p)/(b-a)=0.042 at tau_c=0.044 under constrainAtoB [0,1] - so the effective step
            % was 13x and 24x smaller than the one the unconstrained variant converges with. Dividing
            % by the slope puts the effective step back at the base factor, so this variant INHERITS
            % the stability the unconstrained variant demonstrably has: the multiplier
            % |1-factor*dc/dp| is the same number, not a new one to be argued about. The slopes are
            % evaluated at the answer and vary about 5% between the 1.05x start and there, which would
            % only matter for a price whose multiplier already sat near its boundary of 2.
            % Tr needs no bump: the 0.3 in howtoupdate is already 6x the usual 0.05 precisely because
            % it carries the softplus compensation (see above). And cc=4's log slope is Tr=0.194
            % against softplus's 1-exp(-Tr)=0.176, within 10%, so one factor serves both.
            howtoupdateC=howtoupdate;
            howtoupdateC{1,4}=howtoupdate{1,4}/0.077;  % r:     0.005 -> 0.065
            howtoupdateC{3,4}=howtoupdate{3,4}/0.042;  % tau_c: 0.05  -> 1.19
            if aa==2
                heteroagentoptionsC.fminalgo5.howtoupdate=howtoupdateC;
            elseif aa==3
                heteroagentoptionsC.fminalgo9.howtoupdate=howtoupdateC;
            end
        end
        [p_eqm{aa,cc},GEcondns{aa,cc}]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsA, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsC, simoptions, vfoptions);
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

fprintf('\n=== FHorz PType GEptype: parameter-constraint invariance ===\n')
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
