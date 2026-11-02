function output=CoreFHorzTPath_constraints(T,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline)
% Constraints on the prices of a general eqm transition path: transpathoptions.constrainpositive, constrain0to1
% and constrainAtoB (with constrainAtoBlimits). They are not a change of variables (as in the stationary general
% eqm) but a fraction-to-the-boundary rule applied to each update: a (period,price) that the update would take to
% or past its bound moves only transpathoptions.constraintstepfraction of the way from its current value to the
% bound, and every other element is the update exactly as proposed. Ported from CoreInfHorzTPath_constraints. So:
%   (1) constraints that never bind must give a bit-identical path to no constraints
%   (2) one update with a constraint that binds must equal the rule applied by hand to the unconstrained update;
%   (2b) the same with the step fraction set per price
%   (3) refusals by name, audited as a set
%   (4)-(5) Anderson (GEnewprice=2): non-binding bounds bit-identical; a binding bound keeps the path inside it
%   (6)-(7) the same for the Type-I Anderson, and (8)-(9) for Anderson with local search
%   (10)-(12) quasi-Newton (GEnewprice=1), whose rule shortens the whole Newton step by one scalar: non-binding bounds
%   bit-identical; one step against the rule applied by hand (to r AND w, as the whole step shrinks); a binding bound
%   keeps the path inside it
%   (13) quasi-Newton, every Jacobianmethod, with an upper bound on r closer to the guess than epsprice: the finite
%   differences must perturb inward (one-sided, downward), so the Jacobian, and with it the step, differs from the
%   unconstrained one (the DISCRIMINATOR that the inward branch fired, in the quasi-Newton loop for FullJacobian and
%   LudwigPath, in TransitionPath_FHorz_FullJacobiantriangular, TransitionPath_FHorz_LudwigWstationary and
%   TransitionPath_FHorz_FakeNewsJacobian for the others); and the plain FullJacobian loop and the triangular one must
%   still agree with each other exactly
%   Permanent types are checked in CoreFHorzTPathTests' sibling bank, CoreFHorzTPathPTypeTests (CoreFHorzTPathPType_constraints)
% No figures. Uses the nod z model of the doPart(10) and doPart(14) subcodes, and its Cobb-Douglas general eqm.

vfoptions=struct();
simoptions=struct();

n_d=0; d_grid=[];
ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
jequaloneDist=zeros([n_a,n_z],'gpuArray'); jequaloneDist(1,ceil(n_z/2))=1;

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

FnsToEvaluateGE.K=@(aprime,a,z) a;
FnsToEvaluateGE.N=@(aprime,a,z,kappa_j) kappa_j*z;

GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

heteroagentoptionsGE=struct();
[p_eqm,GEcondns]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluateGE,GeneralEqmEqnsGE,Params,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptionsGE,simoptions,vfoptions);
fprintf('constraints: stationary GE r=%2.8f, w=%2.8f, conditions CapitalMarket=%.3e LabourMarket=%.3e \n',p_eqm.r,p_eqm.w,GEcondns.CapitalMarket,GEcondns.LabourMarket)

ParamsGE=Params;
ParamsGE.r=p_eqm.r;
ParamsGE.w=p_eqm.w;
[V_finalGE,Policy_finalGE]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions);
AgentDist_initialGE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_finalGE,n_d,n_a,n_z,N_j,pi_z,ParamsGE,simoptions);

ParamPathGE.sigma=Params.sigma*ones(1,T); % constant: nothing changes, this is a null reform

transpathoptionsS=transpathoptionsbaseline;
transpathoptionsS.GEnewprice=3; % shooting
transpathoptionsS.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; ...
                                           'LabourMarket','w',0,0.1};
transpathoptionsS.verbose=0;

%% (1) Constraints that never bind give exactly the path without them
% The usual bumped guess (both prices up to 5% off, over the first two-thirds of the path), and a fixed 25
% iterations: the claim is that the rule leaves every update it does not cut untouched, which holds iteration by
% iteration, so a converged solve is not needed. The bounds are far from anything the path visits: r in
% (-firmdelta,0.2) and w positive.
Tbump=floor(2*T/3);
PricePathBumped.r=p_eqm.r*ones(1,T);
PricePathBumped.w=p_eqm.w*ones(1,T);
PricePathBumped.r(1:Tbump)=PricePathBumped.r(1:Tbump).*(1+0.01+0.04*linspace(0,1,Tbump));
PricePathBumped.w(1:Tbump)=PricePathBumped.w(1:Tbump).*(1+0.05-0.04*linspace(0,1,Tbump));
transpathoptions1=transpathoptionsS;
transpathoptions1.maxiter=25;
[PricePath1u,GEcondnPath1u]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions1, simoptions, vfoptions);
transpathoptions1c=transpathoptions1;
transpathoptions1c.constrainpositive={'w'};
transpathoptions1c.constrainAtoB={'r'};
transpathoptions1c.constrainAtoBlimits.r=[-Params.firmdelta,0.2];
[PricePath1c,GEcondnPath1c]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions1c, simoptions, vfoptions);
fprintf('constraints: (1) bounds that never bind, 25 shooting iterations, this should be zero, PricePath r: %.3e, w: %.3e \n',max(abs(PricePath1c.r-PricePath1u.r)),max(abs(PricePath1c.w-PricePath1u.w)))
fprintf('constraints: (1) bounds that never bind, 25 shooting iterations, this should be zero, GEcondnPath CapitalMarket: %.3e, LabourMarket: %.3e \n',max(abs(GEcondnPath1c.CapitalMarket-GEcondnPath1u.CapitalMarket)),max(abs(GEcondnPath1c.LabourMarket-GEcondnPath1u.LabourMarket)))

%% (2) One update with a binding bound is the rule applied by hand to the unconstrained update
% A flat guess, r 5% above its equilibrium in every period (including T, which is never updated, so the guess is
% strictly inside whatever bound is set below). One unconstrained update (maxiter=1 returns the guess after one
% update) shows which way r moves; the bound is then put midway between the smallest and the largest
% move, so it cuts the larger moves and leaves the smaller ones alone. One constrained update from the same guess must then equal,
% element by element, r0+theta*(bound-r0) where the unconstrained update reached or passed the bound, and the
% unconstrained update itself everywhere else. w carries a bound that never binds.
PricePathFlat.r=1.05*p_eqm.r*ones(1,T);
PricePathFlat.w=p_eqm.w*ones(1,T);
r0=PricePathFlat.r(1);
transpathoptions2=transpathoptionsS;
transpathoptions2.maxiter=1;
PricePath2u=TransitionPath_Case1_FHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions2, simoptions, vfoptions);
moves=gather(PricePath2u.r(1:T-1))-r0;
transpathoptions2c=transpathoptions2;
transpathoptions2c.constrainpositive={'w'};
transpathoptions2c.constrainAtoB={'r'};
theta=0.5; % the default constraintstepfraction, left unset in transpathoptions2c so the default is what is checked
% The bound goes midway between the smallest and the largest move (in the direction r moves), not halfway along the
% largest one: in this life-cycle model one update moves r by nearly the same amount in every period, so a bound at
% half the largest move cuts every period and nothing checks what is NOT cut (seen 2026-11-02)
if max(-moves)>=max(moves)
    % r moves down: a lower bound midway between the smallest and the largest downward move
    bound=r0-0.5*(max(-moves)+max(min(-moves),0));
    transpathoptions2c.constrainAtoBlimits.r=[bound,1];
    crossed=(gather(PricePath2u.r(1:T-1))<=bound);
else
    % r moves up: an upper bound midway between the smallest and the largest upward move
    bound=r0+0.5*(max(moves)+max(min(moves),0));
    transpathoptions2c.constrainAtoBlimits.r=[-1,bound];
    crossed=(gather(PricePath2u.r(1:T-1))>=bound);
end
PricePath2c=TransitionPath_Case1_FHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions2c, simoptions, vfoptions);
rbyhand=gather(PricePath2u.r(1:T-1));
rbyhand(crossed)=r0+theta*(bound-r0);
fprintf('constraints: (2) the unconstrained update reaches the bound somewhere (else this tests nothing), this should be one: %i \n',any(crossed))
fprintf('constraints: (2) and leaves it alone somewhere (else nothing checks what is NOT cut), this should be one: %i \n',any(~crossed))
fprintf('constraints: (2) one constrained update vs the rule applied by hand to the unconstrained one, this should be zero, r: %.3e \n',max(abs(gather(PricePath2c.r(1:T-1))-rbyhand)))
fprintf('constraints: (2) w has a bound that does not bind, so its update is untouched, this should be zero: %.3e \n',max(abs(gather(PricePath2c.w-PricePath2u.w))))
fprintf('constraints: (2) period T is never updated, this should be zero: %.3e \n',abs(gather(PricePath2c.r(T))-r0))
fprintf('constraints: (2) the constrained update is strictly inside the bound, this should be one: %i \n',all(abs(gather(PricePath2c.r(1:T-1))-bound)>0 & sign(gather(PricePath2c.r(1:T-1))-bound)==sign(r0-bound)))

% (2b) The same update with the step fraction set per price: constraintstepfraction.r=0.25 (w gets one too, 0.9, but
% its bound does not bind, so that one must change nothing). The rule by hand now uses 0.25.
transpathoptions2d=transpathoptions2c;
transpathoptions2d.constraintstepfraction.r=0.25;
transpathoptions2d.constraintstepfraction.w=0.9;
PricePath2d=TransitionPath_Case1_FHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions2d, simoptions, vfoptions);
rbyhand2d=gather(PricePath2u.r(1:T-1));
rbyhand2d(crossed)=r0+0.25*(bound-r0);
fprintf('constraints: (2b) per-price step fraction (r 0.25), one constrained update vs the rule applied by hand, this should be zero, r: %.3e \n',max(abs(gather(PricePath2d.r(1:T-1))-rbyhand2d)))
fprintf('constraints: (2b) per-price step fraction, w (whose bound does not bind) is untouched, this should be zero: %.3e \n',max(abs(gather(PricePath2d.w-PricePath2u.w))))
fprintf('constraints: (2b) per-price 0.25 vs the default 0.5, this should NOT be zero: %.3e \n',max(abs(gather(PricePath2d.r-PricePath2c.r))))

%% (3) Refusals, by name
% Each is a mistake a user can make, and each would otherwise run (or fail somewhere unhelpful). Audited as a set:
% every line is 'should be one'.
err=zeros(1,8);
% A name that is not a price in PricePath
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.constrainpositive={'sigma'};
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(1)=contains(ME.message,'is not a price in PricePath');
end
% The same price under two kinds of constraint
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.constrainpositive={'w'}; transpathoptionsE.constrain0to1={'w'};
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(2)=contains(ME.message,'appears in more than one');
end
% constrainAtoB without its limits
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.constrainAtoB={'r'};
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(3)=contains(ME.message,'you must set transpathoptions.constrainAtoBlimits.r');
end
% Limits the wrong way round
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.constrainAtoB={'r'}; transpathoptionsE.constrainAtoBlimits.r=[0.2,-0.05];
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(4)=contains(ME.message,'must be [A,B] with A<B');
end
% A step fraction of one would put the price on the bound
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.constrainpositive={'w'}; transpathoptionsE.constraintstepfraction=1;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(5)=contains(ME.message,'constraintstepfraction must be strictly between 0 and 1');
end
% An initial guess outside the bound (the bumped r goes up to 5% above its equilibrium, so an upper bound at the
% equilibrium r excludes it)
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.constrainAtoB={'r'}; transpathoptionsE.constrainAtoBlimits.r=[-0.05,p_eqm.r];
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(6)=contains(ME.message,'is not strictly inside its constraint');
end
% A per-price step fraction for a price that is not constrained
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.constrainpositive={'w'}; transpathoptionsE.constraintstepfraction.r=0.25;
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(7)=contains(ME.message,'is not a constrained price');
end
% A step fraction that is neither a scalar nor a struct by name
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.constrainpositive={'w'}; transpathoptionsE.constraintstepfraction=[0.5,0.5];
try
    TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(8)=contains(ME.message,'must be a scalar, or a struct with one field per constrained price');
end
fprintf('constraints: (3) a constrained name that is not a price is refused by name, should be one: %i \n',err(1))
fprintf('constraints: (3) a price under two kinds of constraint is refused by name, should be one: %i \n',err(2))
fprintf('constraints: (3) constrainAtoB without its limits is refused by name, should be one: %i \n',err(3))
fprintf('constraints: (3) limits with A>=B are refused by name, should be one: %i \n',err(4))
fprintf('constraints: (3) constraintstepfraction=1 is refused by name, should be one: %i \n',err(5))
fprintf('constraints: (3) an initial guess outside its constraint is refused by name, should be one: %i \n',err(6))
fprintf('constraints: (3) a per-price step fraction for an unconstrained price is refused by name, should be one: %i \n',err(7))
fprintf('constraints: (3) a step fraction that is neither scalar nor struct is refused by name, should be one: %i \n',err(8))

%% (4) Anderson: constraints that never bind give exactly the path without them
% As (1), with GEnewprice=2. Ten iterations: the default warmup is two plain steps, so eight of them are Anderson
% steps. The rule on the shooting step and the feasibility test on the Anderson step both leave everything inside
% the bounds alone, so the two runs must be bit-identical.
transpathoptionsA=transpathoptionsS;
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=transpathoptionsS.GEnewprice3.howtoupdate;
transpathoptionsA.maxiter=10;
[PricePath4u,GEcondnPath4u]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsA, simoptions, vfoptions);
transpathoptions4c=transpathoptionsA;
transpathoptions4c.constrainpositive={'w'};
transpathoptions4c.constrainAtoB={'r'};
transpathoptions4c.constrainAtoBlimits.r=[-Params.firmdelta,0.2];
[PricePath4c,GEcondnPath4c]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions4c, simoptions, vfoptions);
fprintf('constraints: (4) Anderson, bounds that never bind, 10 iterations, this should be zero, PricePath r: %.3e, w: %.3e \n',max(abs(PricePath4c.r-PricePath4u.r)),max(abs(PricePath4c.w-PricePath4u.w)))
fprintf('constraints: (4) Anderson, bounds that never bind, 10 iterations, this should be zero, GEcondnPath CapitalMarket: %.3e, LabourMarket: %.3e \n',max(abs(GEcondnPath4c.CapitalMarket-GEcondnPath4u.CapitalMarket)),max(abs(GEcondnPath4c.LabourMarket-GEcondnPath4u.LabourMarket)))

%% (5) Anderson: a binding bound keeps the path inside it
% From the flat guess of (2), six unconstrained iterations show where r goes; a bound halfway along its largest move
% then binds. The constrained run must end strictly inside the bound, where the unconstrained one ended outside or
% on it, and the two must differ. No by-hand oracle: once an Anderson step is not taken, the history restarts and
% the two runs part ways.
transpathoptions5=transpathoptionsA;
transpathoptions5.maxiter=6;
PricePath5u=TransitionPath_Case1_FHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions5, simoptions, vfoptions);
moves5=gather(PricePath5u.r(1:T-1))-r0;
transpathoptions5c=transpathoptions5;
transpathoptions5c.constrainpositive={'w'};
transpathoptions5c.constrainAtoB={'r'};
if max(-moves5)>=max(moves5)
    bound5=r0-0.5*max(-moves5);
    transpathoptions5c.constrainAtoBlimits.r=[bound5,1];
    crossed5=any(gather(PricePath5u.r(1:T-1))<=bound5);
else
    bound5=r0+0.5*max(moves5);
    transpathoptions5c.constrainAtoBlimits.r=[-1,bound5];
    crossed5=any(gather(PricePath5u.r(1:T-1))>=bound5);
end
PricePath5c=TransitionPath_Case1_FHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions5c, simoptions, vfoptions);
fprintf('constraints: (5) Anderson, the unconstrained path ends on or past the bound (else this tests nothing), this should be one: %i \n',crossed5)
fprintf('constraints: (5) Anderson, the constrained path ends strictly inside the bound, this should be one: %i \n',all(sign(gather(PricePath5c.r(1:T-1))-bound5)==sign(r0-bound5)))
fprintf('constraints: (5) Anderson, constrained vs unconstrained path, this should NOT be zero: %.3e \n',max(abs(gather(PricePath5c.r-PricePath5u.r))))

%% (6)-(7) The Type-I Anderson: the same two checks as (4) and (5)
transpathoptionsI=transpathoptionsA;
transpathoptionsI.anderson.type='I';
[PricePath6u,GEcondnPath6u]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsI, simoptions, vfoptions);
transpathoptions6c=transpathoptionsI;
transpathoptions6c.constrainpositive={'w'};
transpathoptions6c.constrainAtoB={'r'};
transpathoptions6c.constrainAtoBlimits.r=[-Params.firmdelta,0.2];
[PricePath6c,GEcondnPath6c]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions6c, simoptions, vfoptions);
fprintf('constraints: (6) Type-I Anderson, bounds that never bind, 10 iterations, this should be zero, PricePath r: %.3e, w: %.3e \n',max(abs(PricePath6c.r-PricePath6u.r)),max(abs(PricePath6c.w-PricePath6u.w)))
fprintf('constraints: (6) Type-I Anderson, bounds that never bind, 10 iterations, this should be zero, GEcondnPath CapitalMarket: %.3e, LabourMarket: %.3e \n',max(abs(GEcondnPath6c.CapitalMarket-GEcondnPath6u.CapitalMarket)),max(abs(GEcondnPath6c.LabourMarket-GEcondnPath6u.LabourMarket)))

transpathoptions7=transpathoptionsI;
transpathoptions7.maxiter=6;
PricePath7u=TransitionPath_Case1_FHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions7, simoptions, vfoptions);
moves7=gather(PricePath7u.r(1:T-1))-r0;
transpathoptions7c=transpathoptions7;
transpathoptions7c.constrainpositive={'w'};
transpathoptions7c.constrainAtoB={'r'};
if max(-moves7)>=max(moves7)
    bound7=r0-0.5*max(-moves7);
    transpathoptions7c.constrainAtoBlimits.r=[bound7,1];
    crossed7=any(gather(PricePath7u.r(1:T-1))<=bound7);
else
    bound7=r0+0.5*max(moves7);
    transpathoptions7c.constrainAtoBlimits.r=[-1,bound7];
    crossed7=any(gather(PricePath7u.r(1:T-1))>=bound7);
end
PricePath7c=TransitionPath_Case1_FHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions7c, simoptions, vfoptions);
fprintf('constraints: (7) Type-I Anderson, the unconstrained path ends on or past the bound (else this tests nothing), this should be one: %i \n',crossed7)
fprintf('constraints: (7) Type-I Anderson, the constrained path ends strictly inside the bound, this should be one: %i \n',all(sign(gather(PricePath7c.r(1:T-1))-bound7)==sign(r0-bound7)))
fprintf('constraints: (7) Type-I Anderson, constrained vs unconstrained path, this should NOT be zero: %.3e \n',max(abs(gather(PricePath7c.r-PricePath7u.r))))

%% (8)-(9) Anderson with local search (AndersonAccelerationTPath): the same two checks
vfoptionsLS=vfoptions;
vfoptionsLS.localsearch=1;
vfoptionsLS.nlocalsearch=3;
[PricePath8u,GEcondnPath8u]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsA, simoptions, vfoptionsLS);
[PricePath8c,GEcondnPath8c]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions4c, simoptions, vfoptionsLS);
fprintf('constraints: (8) Anderson with local search, bounds that never bind, 10 iterations, this should be zero, PricePath r: %.3e, w: %.3e \n',max(abs(PricePath8c.r-PricePath8u.r)),max(abs(PricePath8c.w-PricePath8u.w)))
fprintf('constraints: (8) Anderson with local search, bounds that never bind, 10 iterations, this should be zero, GEcondnPath CapitalMarket: %.3e, LabourMarket: %.3e \n',max(abs(GEcondnPath8c.CapitalMarket-GEcondnPath8u.CapitalMarket)),max(abs(GEcondnPath8c.LabourMarket-GEcondnPath8u.LabourMarket)))

PricePath9u=TransitionPath_Case1_FHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions5, simoptions, vfoptionsLS);
moves9=gather(PricePath9u.r(1:T-1))-r0;
transpathoptions9c=transpathoptions5;
transpathoptions9c.constrainpositive={'w'};
transpathoptions9c.constrainAtoB={'r'};
if max(-moves9)>=max(moves9)
    bound9=r0-0.5*max(-moves9);
    transpathoptions9c.constrainAtoBlimits.r=[bound9,1];
    crossed9=any(gather(PricePath9u.r(1:T-1))<=bound9);
else
    bound9=r0+0.5*max(moves9);
    transpathoptions9c.constrainAtoBlimits.r=[-1,bound9];
    crossed9=any(gather(PricePath9u.r(1:T-1))>=bound9);
end
PricePath9c=TransitionPath_Case1_FHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions9c, simoptions, vfoptionsLS);
fprintf('constraints: (9) Anderson with local search, the unconstrained path ends on or past the bound (else this tests nothing), this should be one: %i \n',crossed9)
fprintf('constraints: (9) Anderson with local search, the constrained path ends strictly inside the bound, this should be one: %i \n',all(sign(gather(PricePath9c.r(1:T-1))-bound9)==sign(r0-bound9)))
fprintf('constraints: (9) Anderson with local search, constrained vs unconstrained path, this should NOT be zero: %.3e \n',max(abs(gather(PricePath9c.r-PricePath9u.r))))

%% (10)-(12) Quasi-Newton (GEnewprice=1)
% The rule here shortens the WHOLE Newton step by one scalar, alpha, so the direction is kept: alpha is set by the
% element whose full step would go furthest past its bound, which moves constraintstepfraction (default 0.99 for
% quasi-Newton) of the way to it. Default Jacobianmethod, LudwigPath.
transpathoptionsQ=transpathoptionsS;
transpathoptionsQ.GEnewprice=1;
transpathoptionsQ.maxiter=3;
% (10) bounds that never bind: bit-identical
[PricePath10u,GEcondnPath10u]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsQ, simoptions, vfoptions);
transpathoptions10c=transpathoptionsQ;
transpathoptions10c.constrainpositive={'w'};
transpathoptions10c.constrainAtoB={'r'};
transpathoptions10c.constrainAtoBlimits.r=[-Params.firmdelta,0.2];
[PricePath10c,GEcondnPath10c]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions10c, simoptions, vfoptions);
fprintf('constraints: (10) quasi-Newton, bounds that never bind, 3 iterations, this should be zero, PricePath r: %.3e, w: %.3e \n',max(abs(PricePath10c.r-PricePath10u.r)),max(abs(PricePath10c.w-PricePath10u.w)))
fprintf('constraints: (10) quasi-Newton, bounds that never bind, 3 iterations, this should be zero, GEcondnPath CapitalMarket: %.3e, LabourMarket: %.3e \n',max(abs(GEcondnPath10c.CapitalMarket-GEcondnPath10u.CapitalMarket)),max(abs(GEcondnPath10c.LabourMarket-GEcondnPath10u.LabourMarket)))

% (11) One Newton step against the rule applied by hand. From the flat guess of (2), one unconstrained step gives
% p0+dx; with a bound halfway along r's largest move, the constrained step must be p0+alpha*dx with
% alpha=0.99*0.5=0.495, for r AND for w, whose own bound does not bind but whose step is shortened all the same.
% That is the whole point of the scalar rule, so w is the half of this check that can tell it from an
% element-by-element one. Valid only if both steps were ACCEPTED (a step that makes the general eqm conditions worse
% is undone and the guess returned unchanged), so both are flagged.
transpathoptions11=transpathoptionsQ;
transpathoptions11.maxiter=1;
PricePath11u=TransitionPath_Case1_FHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions11, simoptions, vfoptions);
dr=gather(PricePath11u.r(1:T-1))-r0;
dw=gather(PricePath11u.w(1:T-1))-gather(PricePathFlat.w(1:T-1));
transpathoptions11c=transpathoptions11;
transpathoptions11c.constrainpositive={'w'};
transpathoptions11c.constrainAtoB={'r'};
if max(-dr)>=max(dr)
    bound11=r0-0.5*max(-dr);
    transpathoptions11c.constrainAtoBlimits.r=[bound11,1];
else
    bound11=r0+0.5*max(dr);
    transpathoptions11c.constrainAtoBlimits.r=[-1,bound11];
end
PricePath11c=TransitionPath_Case1_FHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions11c, simoptions, vfoptions);
alpha11=0.99*0.5;
fprintf('constraints: (11) quasi-Newton, the unconstrained step was accepted (else this tests nothing), this should be one: %i \n',max(abs(dr))>0)
fprintf('constraints: (11) quasi-Newton, the constrained step was accepted (else this tests nothing), this should be one: %i \n',max(abs(gather(PricePath11c.r(1:T-1))-r0))>0)
fprintf('constraints: (11) quasi-Newton, one constrained step vs p0+alpha*dx by hand, this should be zero, r: %.3e \n',max(abs(gather(PricePath11c.r(1:T-1))-(r0+alpha11*dr))))
fprintf('constraints: (11) quasi-Newton, one constrained step vs p0+alpha*dx by hand, w (not binding, but shortened with the rest), this should be zero: %.3e \n',max(abs(gather(PricePath11c.w(1:T-1))-(gather(PricePathFlat.w(1:T-1))+alpha11*dw))))

% (12) A few iterations with the bound of (11): the path stays strictly inside it, and differs from the unconstrained
transpathoptions12=transpathoptionsQ;
transpathoptions12.maxiter=5;
PricePath12u=TransitionPath_Case1_FHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions12, simoptions, vfoptions);
transpathoptions12c=transpathoptions11c;
transpathoptions12c.maxiter=5;
PricePath12c=TransitionPath_Case1_FHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions12c, simoptions, vfoptions);
fprintf('constraints: (12) quasi-Newton, the constrained path stays strictly inside the bound, this should be one: %i \n',all(sign(gather(PricePath12c.r(1:T-1))-bound11)==sign(r0-bound11)))
fprintf('constraints: (12) quasi-Newton, constrained vs unconstrained path, this should NOT be zero: %.3e \n',max(abs(gather(PricePath12c.r-PricePath12u.r))))

%% (13) Quasi-Newton, every Jacobianmethod, with an upper bound on r within epsprice of the guess
% From the flat guess of (2), with r capped at r0+epsprice/2: every finite difference in r would now step past the
% bound, so each Jacobian method must perturb r downward instead (a one-sided difference the other way). r is above
% its equilibrium, so the Newton step moves it down and the bound does not cut the step itself (flagged: if it did,
% the 'NOT zero' below could come from the cut rather than from the inward perturbation). So for each method the
% one constrained step differs from the unconstrained one only through the Jacobian, which is the discriminator that
% its inward branch fired. And the two FullJacobian builds, the plain loop and the triangular one (which has its own
% copy of the inward branch), must still give the same step exactly, as they do without constraints (CoreFHorzTPathAlgo
% check (G)). A shorter path, as FullJacobian is (T-1)*2 path solves per build.
T13=20;
PricePathFlat13.r=PricePathFlat.r(1:T13);
PricePathFlat13.w=PricePathFlat.w(1:T13);
ParamPath13.sigma=Params.sigma*ones(1,T13);
epsprice=1e-3; % the default transpathoptions.epsprice, set explicitly below so that the bound is placed against it
methods13={'FullJacobian','FullJacobian','LudwigPath','LudwigStationary','LudwigSSJ'};
reuse13=[0,1,0,0,0];
stepc13=cell(1,length(methods13));
for mm=1:length(methods13)
    transpathoptions13=transpathoptionsS;
    transpathoptions13.GEnewprice=1;
    transpathoptions13.maxiter=1;
    transpathoptions13.epsprice=epsprice;
    transpathoptions13.GEnewprice1.Jacobianmethod=methods13{mm};
    transpathoptions13.GEnewprice1.FullJacobianReuseVpath=reuse13(mm);
    PricePath13u=TransitionPath_Case1_FHorz(PricePathFlat13, ParamPath13, T13, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions13, simoptions, vfoptions);
    transpathoptions13c=transpathoptions13;
    transpathoptions13c.constrainpositive={'w'};
    transpathoptions13c.constrainAtoB={'r'};
    transpathoptions13c.constrainAtoBlimits.r=[-1,r0+0.5*epsprice];
    PricePath13c=TransitionPath_Case1_FHorz(PricePathFlat13, ParamPath13, T13, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions13c, simoptions, vfoptions);
    stepu=[gather(PricePath13u.r(1:T13-1))-r0, gather(PricePath13u.w(1:T13-1))-gather(PricePathFlat13.w(1:T13-1))];
    stepc=[gather(PricePath13c.r(1:T13-1))-r0, gather(PricePath13c.w(1:T13-1))-gather(PricePathFlat13.w(1:T13-1))];
    stepc13{mm}=stepc;
    lab=sprintf('%s (FullJacobianReuseVpath=%i)',methods13{mm},reuse13(mm));
    fprintf('constraints: (13) quasi-Newton %s, both steps accepted (else this tests nothing), this should be one: %i \n',lab,max(abs(stepu))>0 && max(abs(stepc))>0)
    fprintf('constraints: (13) quasi-Newton %s, the unconstrained step leaves r below the bound, so the bound does not cut the step, this should be one: %i \n',lab,all(gather(PricePath13u.r(1:T13-1))<r0+0.5*epsprice))
    fprintf('constraints: (13) quasi-Newton %s, the inward finite difference changes the step, this should NOT be zero: %.3e \n',lab,max(abs(stepc-stepu)))
end
fprintf('constraints: (13) quasi-Newton FullJacobian with the inward finite difference, plain loop vs triangular, this should be zero: %.3e \n',max(abs(stepc13{2}-stepc13{1})))

output=1;

end
