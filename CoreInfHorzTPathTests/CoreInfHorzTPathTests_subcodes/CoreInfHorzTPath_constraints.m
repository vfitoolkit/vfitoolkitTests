function output=CoreInfHorzTPath_constraints(T,n_a_GE,n_z,a_grid_GE,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline)
% Constraints on the prices of a general eqm transition path: transpathoptions.constrainpositive, constrain0to1
% and constrainAtoB (with constrainAtoBlimits). They are not a change of variables (as in the stationary general
% eqm) but a fraction-to-the-boundary rule applied to each update: a (period,price) that the update would take to
% or past its bound moves only transpathoptions.constraintstepfraction of the way from its current value to the
% bound, and every other element is the update exactly as proposed. So:
%   (1) constraints that never bind must give a bit-identical path to no constraints
%   (2) one update with a constraint that binds must equal the rule applied by hand to the unconstrained update
%   (3) refusals by name
% Shooting (GEnewprice=3) only, for now. Anderson and quasi-Newton get the same checks when they get the rule.
% No figures. Uses the nod z model of the doPart(1) subcodes, and its Cobb-Douglas general eqm.

vfoptions=struct();
simoptions=struct();

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
n_d=0; d_grid=[];

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

FnsToEvaluateGE.K=@(aprime,a,z) a;
FnsToEvaluateGE.N=@(aprime,a,z) z;

GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

heteroagentoptionsGE=struct();
[p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a_GE, n_z, 0, pi_z, d_grid, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptions, vfoptions);
fprintf('constraints: stationary GE r=%2.8f, w=%2.8f, conditions CapitalMarket=%.3e LabourMarket=%.3e \n',p_eqm.r,p_eqm.w,GEcondns.CapitalMarket,GEcondns.LabourMarket)

ParamsGE=Params;
ParamsGE.r=p_eqm.r;
ParamsGE.w=p_eqm.w;
[V_finalGE,Policy_finalGE]=ValueFnIter_InfHorz(n_d,n_a_GE,n_z,d_grid,a_grid_GE,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions);
AgentDist_initialGE=StationaryDist_InfHorz(Policy_finalGE,n_d,n_a_GE,n_z,pi_z,simoptions,ParamsGE,[]);

ParamPathGE.sigma=Params.sigma*ones(1,T); % constant: nothing changes, this is a null reform

transpathoptionsS=transpathoptionsbaseline;
transpathoptionsS.GEnewprice=3; % shooting
transpathoptionsS.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; ...
                                           'LabourMarket','w',0,0.1};
transpathoptionsS.verbose=0;

%% (1) Constraints that never bind give exactly the path without them
% The usual bumped guess (both prices up to 5% off, over the first two-thirds of the path), and a fixed 25
% iterations: the claim is that the rule leaves every update it does not cut untouched, which holds iteration by
% iteration, so a converged solve is not needed (and 25 is cheap). The bounds are far from anything the path visits:
% r in (-firmdelta,0.2) and w positive. NOT r below 1/beta-1, which looks natural but is only 0.0526 here against
% an equilibrium r of about 0.0525, so the bumped guess itself would be outside it.
Tbump=floor(2*T/3);
PricePathBumped.r=p_eqm.r*ones(1,T);
PricePathBumped.w=p_eqm.w*ones(1,T);
PricePathBumped.r(1:Tbump)=PricePathBumped.r(1:Tbump).*(1+0.01+0.04*linspace(0,1,Tbump));
PricePathBumped.w(1:Tbump)=PricePathBumped.w(1:Tbump).*(1+0.05-0.04*linspace(0,1,Tbump));
transpathoptions1=transpathoptionsS;
transpathoptions1.maxiter=25;
[PricePath1u,GEcondnPath1u]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptions1, simoptions, vfoptions);
transpathoptions1c=transpathoptions1;
transpathoptions1c.constrainpositive={'w'};
transpathoptions1c.constrainAtoB={'r'};
transpathoptions1c.constrainAtoBlimits.r=[-Params.firmdelta,0.2];
[PricePath1c,GEcondnPath1c]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptions1c, simoptions, vfoptions);
fprintf('constraints: (1) bounds that never bind, 25 shooting iterations, this should be zero, PricePath r: %.3e, w: %.3e \n',max(abs(PricePath1c.r-PricePath1u.r)),max(abs(PricePath1c.w-PricePath1u.w)))
fprintf('constraints: (1) bounds that never bind, 25 shooting iterations, this should be zero, GEcondnPath CapitalMarket: %.3e, LabourMarket: %.3e \n',max(abs(GEcondnPath1c.CapitalMarket-GEcondnPath1u.CapitalMarket)),max(abs(GEcondnPath1c.LabourMarket-GEcondnPath1u.LabourMarket)))

%% (2) One update with a binding bound is the rule applied by hand to the unconstrained update
% A flat guess, r 5% above its equilibrium in every period (including T, which is never updated, so the guess is
% strictly inside whatever bound is set below). One unconstrained update (maxiter=1 returns the guess after one
% update) shows which way r moves; the bound is then put halfway along the largest move, so it cuts that element
% and possibly others, and leaves the rest alone. One constrained update from the same guess must then equal,
% element by element, r0+theta*(bound-r0) where the unconstrained update reached or passed the bound, and the
% unconstrained update itself everywhere else. That is an exact oracle for the rule, and it checks both what is cut
% and what is not. w carries a bound that never binds, so a column without a binding bound is in there too.
PricePathFlat.r=1.05*p_eqm.r*ones(1,T);
PricePathFlat.w=p_eqm.w*ones(1,T);
r0=PricePathFlat.r(1);
transpathoptions2=transpathoptionsS;
transpathoptions2.maxiter=1;
PricePath2u=TransitionPath_InfHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptions2, simoptions, vfoptions);
moves=gather(PricePath2u.r(1:T-1))-r0;
transpathoptions2c=transpathoptions2;
transpathoptions2c.constrainpositive={'w'};
transpathoptions2c.constrainAtoB={'r'};
theta=0.5; % the default constraintstepfraction, left unset in transpathoptions2c so the default is what is checked
if max(-moves)>=max(moves)
    % r moves down: a lower bound halfway along the largest downward move
    bound=r0-0.5*max(-moves);
    transpathoptions2c.constrainAtoBlimits.r=[bound,1];
    crossed=(gather(PricePath2u.r(1:T-1))<=bound);
else
    % r moves up: an upper bound halfway along the largest upward move
    bound=r0+0.5*max(moves);
    transpathoptions2c.constrainAtoBlimits.r=[-1,bound];
    crossed=(gather(PricePath2u.r(1:T-1))>=bound);
end
PricePath2c=TransitionPath_InfHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptions2c, simoptions, vfoptions);
rbyhand=gather(PricePath2u.r(1:T-1));
rbyhand(crossed)=r0+theta*(bound-r0);
fprintf('constraints: (2) the unconstrained update reaches the bound somewhere (else this tests nothing), this should be one: %i \n',any(crossed))
fprintf('constraints: (2) and leaves it alone somewhere (else nothing checks what is NOT cut), this should be one: %i \n',any(~crossed))
fprintf('constraints: (2) one constrained update vs the rule applied by hand to the unconstrained one, this should be zero, r: %.3e \n',max(abs(gather(PricePath2c.r(1:T-1))-rbyhand)))
fprintf('constraints: (2) w has a bound that does not bind, so its update is untouched, this should be zero: %.3e \n',max(abs(gather(PricePath2c.w-PricePath2u.w))))
fprintf('constraints: (2) period T is never updated, this should be zero: %.3e \n',abs(gather(PricePath2c.r(T))-r0))
fprintf('constraints: (2) the constrained update is strictly inside the bound, this should be one: %i \n',all(abs(gather(PricePath2c.r(1:T-1))-bound)>0 & sign(gather(PricePath2c.r(1:T-1))-bound)==sign(r0-bound)))

%% (3) Refusals, by name
% Each is a mistake a user can make, and each would otherwise run (or fail somewhere unhelpful)
err=zeros(1,7);
% A name that is not a price in PricePath
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.constrainpositive={'sigma'};
try
    TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(1)=contains(ME.message,'is not a price in PricePath');
end
% The same price under two kinds of constraint
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.constrainpositive={'w'}; transpathoptionsE.constrain0to1={'w'};
try
    TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(2)=contains(ME.message,'appears in more than one');
end
% constrainAtoB without its limits
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.constrainAtoB={'r'};
try
    TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(3)=contains(ME.message,'you must set transpathoptions.constrainAtoBlimits.r');
end
% Limits the wrong way round
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.constrainAtoB={'r'}; transpathoptionsE.constrainAtoBlimits.r=[0.2,-0.05];
try
    TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(4)=contains(ME.message,'must be [A,B] with A<B');
end
% A step fraction of one would put the price on the bound
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.constrainpositive={'w'}; transpathoptionsE.constraintstepfraction=1;
try
    TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(5)=contains(ME.message,'constraintstepfraction must be strictly between 0 and 1');
end
% An initial guess outside the bound (the bumped r goes up to 5% above its equilibrium, so an upper bound at the
% equilibrium r excludes it)
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.constrainAtoB={'r'}; transpathoptionsE.constrainAtoBlimits.r=[-0.05,p_eqm.r];
try
    TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(6)=contains(ME.message,'is not strictly inside its constraint');
end
% Anderson does not apply the rule yet, so it must refuse rather than ignore the constraints. THIS CHECK GOES
% when Anderson gets the rule (stage C2), and the same for quasi-Newton (C3).
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.GEnewprice=2; transpathoptionsE.GEnewprice2.howtoupdate=transpathoptionsS.GEnewprice3.howtoupdate;
transpathoptionsE.constrainpositive={'w'};
try
    TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(7)=contains(ME.message,'so far only implemented for transpathoptions.GEnewprice=3');
end
fprintf('constraints: (3) a constrained name that is not a price is refused by name, should be one: %i \n',err(1))
fprintf('constraints: (3) a price under two kinds of constraint is refused by name, should be one: %i \n',err(2))
fprintf('constraints: (3) constrainAtoB without its limits is refused by name, should be one: %i \n',err(3))
fprintf('constraints: (3) limits with A>=B are refused by name, should be one: %i \n',err(4))
fprintf('constraints: (3) constraintstepfraction=1 is refused by name, should be one: %i \n',err(5))
fprintf('constraints: (3) an initial guess outside its constraint is refused by name, should be one: %i \n',err(6))
fprintf('constraints: (3) constraints with Anderson (not yet implemented) are refused by name, should be one: %i \n',err(7))

output=1;

end
