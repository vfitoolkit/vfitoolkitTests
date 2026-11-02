function output=CoreFHorzTPathPType_constraints(T,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline)
% Constraints on the prices of a general eqm transition path with permanent types (TransitionPath_Case1_FHorz_PType):
% transpathoptions.constrainpositive, constrain0to1 and constrainAtoB (with constrainAtoBlimits), the
% fraction-to-the-boundary rule applied to each update. The version without permanent types, with every algorithm, is
% CoreFHorzTPath_constraints (CoreFHorzTPathTests doPart(15)); this checks that the PType top level sets the bounds and
% that its shooting and Anderson apply them. Moved here from CoreFHorzTPath_constraints 2026-11-03.
%   (1) one shooting update with a binding bound against the rule applied by hand to the unconstrained update, and w
%       (whose bound does not bind) untouched
%   (2) Anderson with bounds that never bind gives a bit-identical path
%   (3) a constrained name that is not a price is refused by name
% Two identical types (ptweights 0.6 and 0.4), each the nod z model of CoreFHorzTPath_constraints, with its
% Cobb-Douglas general eqm in r and w (common to both types). No figures.

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

% The stationary general eqm of one type (the types are identical, so it is that of both)
heteroagentoptionsGE=struct();
[p_eqm,GEcondns]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluateGE,GeneralEqmEqnsGE,Params,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptionsGE,simoptions,vfoptions);
fprintf('PType constraints: stationary GE r=%2.8f, w=%2.8f, conditions CapitalMarket=%.3e LabourMarket=%.3e \n',p_eqm.r,p_eqm.w,GEcondns.CapitalMarket,GEcondns.LabourMarket)

ParamsGE=Params;
ParamsGE.r=p_eqm.r;
ParamsGE.w=p_eqm.w;
[V_finalGE,Policy_finalGE]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions);
AgentDist_initialGE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_finalGE,n_d,n_a,n_z,N_j,pi_z,ParamsGE,simoptions);

ParamPathGE.sigma=Params.sigma*ones(1,T); % constant: nothing changes, this is a null reform

% Two identical types
PTypeDistParamNames={'ptweights_con'};
ParamsPT=ParamsGE;
ParamsPT.ptweights_con=[0.6;0.4];
V_finalPT.ptype001=V_finalGE;
V_finalPT.ptype002=V_finalGE;
AgentDist_initialPT.ptype001=AgentDist_initialGE;
AgentDist_initialPT.ptype002=AgentDist_initialGE;
AgentDist_initialPT.ptweights=[0.6;0.4];

transpathoptionsS=transpathoptionsbaseline;
transpathoptionsS.GEnewprice=3; % shooting
transpathoptionsS.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; ...
                                           'LabourMarket','w',0,0.1};
transpathoptionsS.verbose=0;

%% (1) One shooting update with a binding bound is the rule applied by hand to the unconstrained update
% A flat guess, r 5% above its equilibrium in every period (including T, which is never updated, so the guess is
% strictly inside whatever bound is set below). One unconstrained update (maxiter=1 returns the guess after one
% update) shows which way r moves. The bound goes midway between the smallest and the largest move in that direction
% (not halfway along the largest one: in this life-cycle model one update moves r by nearly the same amount in every
% period, so that would cut every period and nothing would check what is NOT cut). One constrained update from the same
% guess must then equal, element by element, r0+theta*(bound-r0) where the unconstrained update reached or passed the
% bound, and the unconstrained update itself everywhere else.
PricePathFlat.r=1.05*p_eqm.r*ones(1,T);
PricePathFlat.w=p_eqm.w*ones(1,T);
r0=PricePathFlat.r(1);
theta=0.5; % the default constraintstepfraction, left unset so the default is what is checked
transpathoptions1=transpathoptionsS;
transpathoptions1.maxiter=1;
PricePath1u=TransitionPath_Case1_FHorz_PType(PricePathFlat, ParamPathGE, T, V_finalPT, AgentDist_initialPT, jequaloneDist, n_d, n_a, n_z, N_j, 2, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsPT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions1, simoptions, vfoptions);
moves=gather(PricePath1u.r(:)')-r0;
moves=moves(1:T-1);
transpathoptions1c=transpathoptions1;
transpathoptions1c.constrainpositive={'w'};
transpathoptions1c.constrainAtoB={'r'};
if max(-moves)>=max(moves)
    % r moves down: a lower bound midway between the smallest and the largest downward move
    bound=r0-0.5*(max(-moves)+max(min(-moves),0));
    transpathoptions1c.constrainAtoBlimits.r=[bound,1];
    crossed=(r0+moves<=bound);
else
    % r moves up: an upper bound midway between the smallest and the largest upward move
    bound=r0+0.5*(max(moves)+max(min(moves),0));
    transpathoptions1c.constrainAtoBlimits.r=[-1,bound];
    crossed=(r0+moves>=bound);
end
PricePath1c=TransitionPath_Case1_FHorz_PType(PricePathFlat, ParamPathGE, T, V_finalPT, AgentDist_initialPT, jequaloneDist, n_d, n_a, n_z, N_j, 2, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsPT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions1c, simoptions, vfoptions);
rbyhand=r0+moves;
rbyhand(crossed)=r0+theta*(bound-r0);
r1c=gather(PricePath1c.r(:)');
fprintf('PType constraints: (1) the unconstrained update reaches the bound somewhere (else this tests nothing), this should be one: %i \n',any(crossed))
fprintf('PType constraints: (1) and leaves it alone somewhere (else nothing checks what is NOT cut), this should be one: %i \n',any(~crossed))
fprintf('PType constraints: (1) one constrained shooting update vs the rule applied by hand, this should be zero, r: %.3e \n',max(abs(r1c(1:T-1)-rbyhand)))
fprintf('PType constraints: (1) w (whose bound does not bind) is untouched, this should be zero: %.3e \n',max(abs(gather(PricePath1c.w(:)-PricePath1u.w(:)))))

%% (2) Anderson with bounds that never bind is bit-identical
% The usual bumped guess (both prices up to 5% off, over the first two-thirds of the path), ten iterations (two plain
% warmup steps, then Anderson steps); bounds far from anything the path visits: r in (-firmdelta,0.2), w positive.
Tbump=floor(2*T/3);
PricePathBumped.r=p_eqm.r*ones(1,T);
PricePathBumped.w=p_eqm.w*ones(1,T);
PricePathBumped.r(1:Tbump)=PricePathBumped.r(1:Tbump).*(1+0.01+0.04*linspace(0,1,Tbump));
PricePathBumped.w(1:Tbump)=PricePathBumped.w(1:Tbump).*(1+0.05-0.04*linspace(0,1,Tbump));
transpathoptionsA=transpathoptionsS;
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=transpathoptionsS.GEnewprice3.howtoupdate;
transpathoptionsA.maxiter=10;
transpathoptionsAc=transpathoptionsA;
transpathoptionsAc.constrainpositive={'w'};
transpathoptionsAc.constrainAtoB={'r'};
transpathoptionsAc.constrainAtoBlimits.r=[-Params.firmdelta,0.2];
[PricePath2u,GEcondnPath2u]=TransitionPath_Case1_FHorz_PType(PricePathBumped, ParamPathGE, T, V_finalPT, AgentDist_initialPT, jequaloneDist, n_d, n_a, n_z, N_j, 2, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsPT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsA, simoptions, vfoptions);
[PricePath2c,GEcondnPath2c]=TransitionPath_Case1_FHorz_PType(PricePathBumped, ParamPathGE, T, V_finalPT, AgentDist_initialPT, jequaloneDist, n_d, n_a, n_z, N_j, 2, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsPT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsAc, simoptions, vfoptions);
fprintf('PType constraints: (2) Anderson, bounds that never bind, 10 iterations, this should be zero, PricePath r: %.3e, w: %.3e \n',max(abs(gather(PricePath2c.r(:)-PricePath2u.r(:)))),max(abs(gather(PricePath2c.w(:)-PricePath2u.w(:)))))
fprintf('PType constraints: (2) Anderson, bounds that never bind, 10 iterations, this should be zero, GEcondnPath CapitalMarket: %.3e, LabourMarket: %.3e \n',max(abs(gather(GEcondnPath2c.CapitalMarket(:)-GEcondnPath2u.CapitalMarket(:)))),max(abs(gather(GEcondnPath2c.LabourMarket(:)-GEcondnPath2u.LabourMarket(:)))))

%% (3) A constrained name that is not a price is refused by name
errPT=0;
transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
transpathoptionsE.constrainpositive={'sigma'};
try
    TransitionPath_Case1_FHorz_PType(PricePathBumped, ParamPathGE, T, V_finalPT, AgentDist_initialPT, jequaloneDist, n_d, n_a, n_z, N_j, 2, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsPT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    errPT=contains(ME.message,'is not a price in PricePath');
end
fprintf('PType constraints: (3) a constrained name that is not a price is refused by name, should be one: %i \n',errPT)

output=1;

end
