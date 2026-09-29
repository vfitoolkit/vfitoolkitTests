function output=CoreStationaryGE_FHorz_PType_GEptype(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% FHorz stationary GE with permanent types: general eqm conditions that hold CONDITIONAL ON THE
% PERMANENT TYPE (heteroagentoptions.GEptype), which nothing in this bank has ever exercised.
%
% WHAT THE MODEL BECOMES. Nothing in the model changes. Tr (the lump-sum transfer) is handed in as a
% struct with one field per permanent type instead of a scalar, which is what makes it a per-type GE
% price, and GovBudget is named in heteroagentoptions.GEptype, which is what makes it a per-type GE
% condition: each type's labor-tax revenue is rebated to that type, tau*w*N_i=Tr_i. The GE eqn TEXT
% does not change at all - for a GEptype condition the toolkit rewrites each input name that is an
% AggVar or a ptype-dependent parameter to its '_ptypeNNN' version, so 'N' and 'Tr' in the eqn
% already mean this type's own N and Tr. The households of type i see their own Tr because a
% struct-valued parameter is peeled per type in the ordinary way.
%
% Names_i is given here as a CELL of names rather than as the number of types (which is what the
% other parts of this bank use), because everything here is indexed by type: the per-type prices
% going in, and p_eqm.Tr.ptypeNNN coming back out.
%
% THE FIVE CHECKS.
%   G1  Two IDENTICAL types, general eqm conditions evaluated at the initial prices only
%       (heteroagentoptions.maxiter=0, one model solve, no search). With identical types and equal
%       masses N_ptype001=N_ptype002=N, so each per-type GovBudget is arithmetically the
%       economy-wide one: the GEptype run must reproduce the economy-wide run BIT-IDENTICALLY.
%   G2  The same identical-types model, solved. The two Tr must come out equal to each other, and
%       equal to the economy-wide Tr of the ordinary PType solve.
%   G3  Types that actually differ (sigma 2.2 and 1.8), so the two transfers genuinely differ. The
%       per-type conditions must be satisfied at the answer; the two Tr must NOT be equal (without
%       that, every check here could pass while the feature did nothing); and the per-type budget is
%       verified INDEPENDENTLY, by recomputing each type's N from its own agent distribution at the
%       returned prices rather than trusting the solver's own evaluation of its conditions.
%   G4  GEptype on ALL THREE conditions, so r, Tr and tau_c are each per-type. The two types are
%       then two separate economies that happen to be solved together, and the answer must equal two
%       independent one-type (non-PType) GE solves. Note the GE conditions make this exact: the
%       capital market condition is scale free, and the two budget conditions compare a per-type
%       per-capita AggVar against the same parameters as the solo economy does.
%       The GEptype solve is started from the solo answers, which is what keeps a six-unknown
%       fminsearch affordable; landing and staying there is the fixed-point property being tested.
%   G5  A per-type INTERMEDIATE eqn (heteroagentoptions.intermediateEqnsptype): the tax base is
%       computed per type as an intermediate and GovBudget then just compares it to Tr. Same model
%       as G3 written a second way, so at the initial prices the conditions must be bit-identical.
%
% GE prices: r (capital market), Tr (gov budget), tau_c (consumption-tax budget).
% w is hardcoded from r via the firm FOC, so it is NOT a GE price.

n_p=0;
vfoptions.divideandconquer=1; % finite-horizon tests use divide-and-conquer

% Permanent types: two types differing in sigma
N_i=2;
Names_i={'ptype001','ptype002'}; % named, because everything here is indexed by type
sigmabyptype=[2.2,1.8];
Params.sigma=sigmabyptype;    % differs by permanent type
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

heteroagentoptions0=heteroagentoptions; % just evaluate the GE conditions at the prices in Params
heteroagentoptions0.fminalgo=1;
heteroagentoptions0.maxiter=0;
heteroagentoptionsS=heteroagentoptions; % a real solve
heteroagentoptionsS.fminalgo=1;

% Tr as a per-type price, at the same initial guess for both types
Trbyptype.ptype001=Params.Tr;
Trbyptype.ptype002=Params.Tr;

%% G1: identical types, GE conditions at the initial prices, economy-wide vs by-ptype GovBudget
ParamsID=Params;
ParamsID.sigma=2; % identical types
[~,GEcondnsG1base]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsID, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptions, vfoptions);

ParamsG1=ParamsID;
ParamsG1.Tr=Trbyptype;
heteroagentoptionsG1=heteroagentoptions0;
heteroagentoptionsG1.GEptype={'GovBudget'};
[~,GEcondnsG1]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsG1, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsG1, simoptions, vfoptions);

%% G2: the same identical-types model, solved both ways
[p_eqmG2plain,~]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsID, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsS, simoptions, vfoptions);

heteroagentoptionsG2=heteroagentoptionsS;
heteroagentoptionsG2.GEptype={'GovBudget'};
[p_eqmG2,~]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsG1, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsG2, simoptions, vfoptions);

%% G3: types that differ, so the two transfers differ too
ParamsG3=Params; % sigma=[2.2,1.8]
ParamsG3.Tr=Trbyptype;
[p_eqmG3,GEcondnsG3]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsG3, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsG2, simoptions, vfoptions);

% Independent verification of the per-type budget: solve the model at the returned prices and
% recompute each type's labor supply from its own agent distribution, rather than reading the
% solver's own evaluation of the conditions it was minimizing.
ParamsG3at=ParamsG3;
ParamsG3at.r=p_eqmG3.r;
ParamsG3at.tau_c=p_eqmG3.tau_c;
ParamsG3at.Tr=p_eqmG3.Tr;
[~,PolicyG3]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsG3at,DiscountFactorParamNames,vfoptions);
StationaryDistG3=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,PolicyG3,n_d,n_a,n_z,N_j,Names_i,pi_z,ParamsG3at,simoptions);
AggVarsG3=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDistG3, PolicyG3, FnsToEvaluate, ParamsG3at,n_d,n_a,n_z,N_j,Names_i,d_grid, a_grid, z_grid, simoptions);
wG3=(1-Params.alpha)*Params.A*((p_eqmG3.r+Params.delta)/(Params.alpha*Params.A))^(Params.alpha/(Params.alpha-1));
budgetresid=zeros(1,N_i);
for ii=1:N_i
    budgetresid(ii)=Params.tau*wG3*AggVarsG3.N.(Names_i{ii}).Mean-p_eqmG3.Tr.(Names_i{ii});
end

%% G4: GEptype on all three conditions, so the two types are two separate economies
% First the reference: two independent one-type solves, using the non-PType command.
p_eqmsolo=cell(1,N_i);
for ii=1:N_i
    Paramssolo=Params;
    Paramssolo.sigma=sigmabyptype(ii); % one type, so a scalar
    Paramssolo=rmfield(Paramssolo,'ptypemass');
    [p_eqmsolo{ii},~]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Paramssolo, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptionsS, simoptions, vfoptions);
end

% Now the by-ptype solve, started from those answers: six unknowns, and what is being tested is
% that this IS the fixed point, not how far fminsearch can walk to find it.
% Built in temporaries and then assigned whole: Params.r, Params.Tr and Params.tau_c are doubles at
% this point, and a field of a double cannot be dot-indexed into existence.
for ii=1:N_i
    rG4.(Names_i{ii})=p_eqmsolo{ii}.r;
    TrG4.(Names_i{ii})=p_eqmsolo{ii}.Tr;
    tau_cG4.(Names_i{ii})=p_eqmsolo{ii}.tau_c;
end
ParamsG4=Params;
ParamsG4.r=rG4;
ParamsG4.Tr=TrG4;
ParamsG4.tau_c=tau_cG4;
heteroagentoptionsG4=heteroagentoptionsS;
heteroagentoptionsG4.GEptype={'CapitalMarket','GovBudget','ConsTax'};
[p_eqmG4,GEcondnsG4]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsG4, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsG4, simoptions, vfoptions);

%% G5: the per-type tax base as a per-type intermediate eqn
% Same model as G3, written a second way: an intermediate eqn marked as conditional on ptype gets
% its AggVar inputs rewritten to this type's own, and its own name is then available to the GE eqn
% as a per-type quantity too. At the initial prices this must match G3's setup exactly.
[~,GEcondnsG5inline]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsG3, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsG1, simoptions, vfoptions);

GeneralEqmEqnsG5=GeneralEqmEqns;
GeneralEqmEqnsG5.GovBudget=@(taxbase,Tr) taxbase-Tr;
heteroagentoptionsG5=heteroagentoptionsG1; % maxiter=0, GEptype on GovBudget
heteroagentoptionsG5.intermediateEqns.taxbase=@(tau,r,N,alpha,delta,A) tau*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*N;
heteroagentoptionsG5.intermediateEqnsptype={'taxbase'};
[~,GEcondnsG5]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqnsG5, ParamsG3, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsG5, simoptions, vfoptions);

%% Compare
% G1: the per-type GovBudget entries must each equal the economy-wide one, and the other two
% conditions must be untouched.
dG1=max(abs([GEcondnsG1base.CapitalMarket-GEcondnsG1.CapitalMarket, ...
             GEcondnsG1base.ConsTax-GEcondnsG1.ConsTax, ...
             GEcondnsG1base.GovBudget-GEcondnsG1.GovBudget(1), ...
             GEcondnsG1base.GovBudget-GEcondnsG1.GovBudget(2)]));
% G2: identical types, so the two transfers coincide and reproduce the economy-wide solve
dG2a=abs(p_eqmG2.Tr.(Names_i{1})-p_eqmG2.Tr.(Names_i{2}));
dG2b=max(abs([p_eqmG2.r-p_eqmG2plain.r, p_eqmG2.tau_c-p_eqmG2plain.tau_c, ...
              p_eqmG2.Tr.(Names_i{1})-p_eqmG2plain.Tr, p_eqmG2.Tr.(Names_i{2})-p_eqmG2plain.Tr]));
% G3: conditions satisfied, transfers genuinely different, budget verified independently
dG3a=max(abs([GEcondnsG3.CapitalMarket,GEcondnsG3.ConsTax,GEcondnsG3.GovBudget]));
dG3b=abs(p_eqmG3.Tr.(Names_i{1})-p_eqmG3.Tr.(Names_i{2}));
dG3c=max(abs(budgetresid));
% G4: two separate economies, so the by-ptype answer is the two solo answers
dG4a=0;
for ii=1:N_i
    dG4a=max([dG4a, abs(p_eqmG4.r.(Names_i{ii})-p_eqmsolo{ii}.r), ...
                    abs(p_eqmG4.Tr.(Names_i{ii})-p_eqmsolo{ii}.Tr), ...
                    abs(p_eqmG4.tau_c.(Names_i{ii})-p_eqmsolo{ii}.tau_c)]);
end
dG4b=max(abs([p_eqmsolo{1}.r-p_eqmsolo{2}.r, p_eqmsolo{1}.Tr-p_eqmsolo{2}.Tr, p_eqmsolo{1}.tau_c-p_eqmsolo{2}.tau_c]));
% G5: the per-type intermediate against the same thing written inline
dG5=max(abs([GEcondnsG5inline.CapitalMarket-GEcondnsG5.CapitalMarket, ...
             GEcondnsG5inline.ConsTax-GEcondnsG5.ConsTax, ...
             GEcondnsG5inline.GovBudget-GEcondnsG5.GovBudget]));

fprintf('\n=== FHorz PType: general eqm conditions conditional on permanent type (GEptype) ===\n')
fprintf('identical types, by-ptype GovBudget vs economy-wide, at the initial prices, this should be zero: %.3e \n',dG1)
fprintf('identical types solved: the two transfers, this should be near zero: %.3e \n',dG2a)
fprintf('identical types solved: by-ptype vs economy-wide prices, this should be near zero: %.3e \n',dG2b)
fprintf('types differ: r=%.6f tau_c=%.6f Tr=(%.6f, %.6f) \n',p_eqmG3.r,p_eqmG3.tau_c,p_eqmG3.Tr.(Names_i{1}),p_eqmG3.Tr.(Names_i{2}))
fprintf('types differ: general eqm conditions at the answer, this should be near zero: %.3e \n',dG3a)
fprintf('types differ: the two transfers, this should NOT be zero: %.3e \n',dG3b)
fprintf('types differ: tau*w*N_i-Tr_i recomputed from each type own agent dist, this should be near zero: %.3e \n',dG3c)
fprintf('all three conditions by ptype: solo answers r=(%.6f, %.6f) Tr=(%.6f, %.6f) tau_c=(%.6f, %.6f) \n',p_eqmsolo{1}.r,p_eqmsolo{2}.r,p_eqmsolo{1}.Tr,p_eqmsolo{2}.Tr,p_eqmsolo{1}.tau_c,p_eqmsolo{2}.tau_c)
fprintf('all three conditions by ptype vs two solo one-type solves, this should be near zero: %.3e \n',dG4a)
fprintf('the two solo one-type solves, this should NOT be zero: %.3e \n',dG4b)
fprintf('per-type intermediate eqn vs the same tax base inline, this should be zero: %.3e \n',dG5)

output.GEcondnsG1base=GEcondnsG1base;
output.GEcondnsG1=GEcondnsG1;
output.p_eqmG2plain=p_eqmG2plain;
output.p_eqmG2=p_eqmG2;
output.p_eqmG3=p_eqmG3;
output.GEcondnsG3=GEcondnsG3;
output.AggVarsG3=AggVarsG3;
output.budgetresid=budgetresid;
output.p_eqmsolo=p_eqmsolo;
output.p_eqmG4=p_eqmG4;
output.GEcondnsG4=GEcondnsG4;
output.GEcondnsG5inline=GEcondnsG5inline;
output.GEcondnsG5=GEcondnsG5;

end
