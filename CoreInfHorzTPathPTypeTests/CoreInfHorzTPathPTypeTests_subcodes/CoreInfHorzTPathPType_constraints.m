function output=CoreInfHorzTPathPType_constraints(T,n_a_GE,n_z,a_grid_GE,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline)
% Constraints on the prices of the general eqm transition path with permanent types: transpathoptions.constrainpositive,
% constrain0to1 and constrainAtoB (with constrainAtoBlimits), and constraintstepfraction. The rule is that of
% TransitionPath_InfHorz (tested by CoreInfHorzTPath_constraints): any (period,price) that an update would take to or
% past its bound moves only constraintstepfraction of the way from its current value to the bound. A price that depends
% on ptype has the same bounds for every ptype. Every check has an exact oracle:
%   (1) N_i=1 and N_i=2 identical types vs the solver without permanent types, with the same constraints: bounds that
%       never bind (25 shooting iterations); a bound that binds, one shooting update and 25 iterations; and a per-price
%       step fraction
%   (2) the same for Anderson (GEnewprice=2), Type-II and Type-I, 10 iterations, against a binding bound
%   (3) a price by ptype under a binding bound: two 'islands' (each ptype with its own firm, r and w, general eqm by
%       ptype) against two solo constrained solves, by shooting; Anderson with memory=0 against shooting; and Anderson
%       with its default memory inside the bound and different from Anderson without it
%   (4) refusals by name, including an initial guess that is outside its bound for only one ptype
% Every binding case prints that the unconstrained update did reach the bound (else it tests nothing), and that the
% constrained path differs from the unconstrained one ('should NOT be zero').
% No figures. The model is that of CoreInfHorzTPath_constraints: nod z, Cobb-Douglas general eqm on r and w.

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

FnsToEvaluateGE.K=@(aprime,a,z) a;
FnsToEvaluateGE.N=@(aprime,a,z) z;

GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

heteroagentoptionsGE=struct();
[p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a_GE, n_z, 0, pi_z, d_grid, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptions, vfoptions);
fprintf('PType constraints: stationary GE r=%2.8f, w=%2.8f, conditions CapitalMarket=%.3e LabourMarket=%.3e \n',p_eqm.r,p_eqm.w,GEcondns.CapitalMarket,GEcondns.LabourMarket)

ParamsGE=Params;
ParamsGE.r=p_eqm.r;
ParamsGE.w=p_eqm.w;
[V_finalGE,Policy_finalGE]=ValueFnIter_InfHorz(n_d,n_a_GE,n_z,d_grid,a_grid_GE,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions);
AgentDist_initialGE=StationaryDist_InfHorz(Policy_finalGE,n_d,n_a_GE,n_z,pi_z,simoptions,ParamsGE,[]);

ParamPathGE.sigma=Params.sigma*ones(1,T); % constant: nothing changes, this is a null reform

% N_i=1, and N_i=2 identical types (ptweights 0.6 and 0.4)
PTypeDistParamNames={'ptweights_con'};
V_final1.ptype001=V_finalGE;
AgentDist_initial1.ptype001=AgentDist_initialGE;
AgentDist_initial1.ptweights=1;
Params1=ParamsGE;
Params1.ptweights_con=1;
V_final2.ptype001=V_finalGE; V_final2.ptype002=V_finalGE;
AgentDist_initial2.ptype001=AgentDist_initialGE; AgentDist_initial2.ptype002=AgentDist_initialGE;
AgentDist_initial2.ptweights=[0.6;0.4];
Params2=ParamsGE;
Params2.ptweights_con=[0.6;0.4];

transpathoptionsS=transpathoptionsbaseline;
transpathoptionsS.GEnewprice=3; % shooting
transpathoptionsS.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; ...
                                           'LabourMarket','w',0,0.1};
transpathoptionsS.verbose=0;

% The guesses of CoreInfHorzTPath_constraints: bumped (both prices up to 5% off, over the first two-thirds of the
% path), and flat (r 5% above its equilibrium in every period, including T, which is never updated)
Tbump=floor(2*T/3);
PricePathBumped.r=p_eqm.r*ones(1,T);
PricePathBumped.w=p_eqm.w*ones(1,T);
PricePathBumped.r(1:Tbump)=PricePathBumped.r(1:Tbump).*(1+0.01+0.04*linspace(0,1,Tbump));
PricePathBumped.w(1:Tbump)=PricePathBumped.w(1:Tbump).*(1+0.05-0.04*linspace(0,1,Tbump));
PricePathFlat.r=1.05*p_eqm.r*ones(1,T);
PricePathFlat.w=p_eqm.w*ones(1,T);
r0=PricePathFlat.r(1);

%% (1) N_i=1 and N_i=2 identical types vs no-PType, shooting
% (1a) Bounds that never bind: the constrained PType path is bit-identical to the unconstrained one, 25 iterations
transpathoptions1=transpathoptionsS;
transpathoptions1.maxiter=25;
transpathoptions1c=transpathoptions1;
transpathoptions1c.constrainpositive={'w'};
transpathoptions1c.constrainAtoB={'r'};
transpathoptions1c.constrainAtoBlimits.r=[-Params.firmdelta,0.2];
[PricePath1u,GEcondnPath1u]=TransitionPath_InfHorz_PType(PricePathBumped, ParamPathGE, T, V_final2, AgentDist_initial2, n_d,n_a_GE,n_z, 2, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params2, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions1, simoptions, vfoptions);
[PricePath1c,GEcondnPath1c]=TransitionPath_InfHorz_PType(PricePathBumped, ParamPathGE, T, V_final2, AgentDist_initial2, n_d,n_a_GE,n_z, 2, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params2, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions1c, simoptions, vfoptions);
fprintf('PType constraints: (1a) N_i=2, bounds that never bind, 25 shooting iterations, this should be zero, PricePath r: %.3e, w: %.3e \n',max(abs(gather(PricePath1c.r-PricePath1u.r))),max(abs(gather(PricePath1c.w-PricePath1u.w))))
fprintf('PType constraints: (1a) N_i=2, bounds that never bind, 25 shooting iterations, this should be zero, GEcondnPath CapitalMarket: %.3e, LabourMarket: %.3e \n',max(abs(gather(GEcondnPath1c.CapitalMarket-GEcondnPath1u.CapitalMarket))),max(abs(gather(GEcondnPath1c.LabourMarket-GEcondnPath1u.LabourMarket))))

% (1b) A bound that binds, one update: the bound is put halfway along the largest move of r in one unconstrained
% update without permanent types (as in CoreInfHorzTPath_constraints (2)), and w has a bound that never binds
transpathoptions2=transpathoptionsS;
transpathoptions2.maxiter=1; % returns the guess after one update
PricePath2u=TransitionPath_InfHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptions2, simoptions, vfoptions, []);
moves=gather(PricePath2u.r(1:T-1))-r0;
transpathoptions2c=transpathoptions2;
transpathoptions2c.constrainpositive={'w'};
transpathoptions2c.constrainAtoB={'r'};
if max(-moves)>=max(moves)
    bound=r0-0.5*max(-moves);
    transpathoptions2c.constrainAtoBlimits.r=[bound,1];
    crossed=(gather(PricePath2u.r(1:T-1))<=bound);
else
    bound=r0+0.5*max(moves);
    transpathoptions2c.constrainAtoBlimits.r=[-1,bound];
    crossed=(gather(PricePath2u.r(1:T-1))>=bound);
end
[PricePath2c_NoPT,GEcondnPath2c_NoPT]=TransitionPath_InfHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptions2c, simoptions, vfoptions, []);
[PricePath2c_PT1,GEcondnPath2c_PT1]=TransitionPath_InfHorz_PType(PricePathFlat, ParamPathGE, T, V_final1, AgentDist_initial1, n_d,n_a_GE,n_z, 1, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params1, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions2c, simoptions, vfoptions);
[PricePath2c_PT2,GEcondnPath2c_PT2]=TransitionPath_InfHorz_PType(PricePathFlat, ParamPathGE, T, V_final2, AgentDist_initial2, n_d,n_a_GE,n_z, 2, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params2, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions2c, simoptions, vfoptions);
rbyhand=gather(PricePath2u.r(1:T-1));
rbyhand(crossed)=r0+0.5*(bound-r0); % 0.5 is the default constraintstepfraction
fprintf('PType constraints: (1b) the unconstrained update reaches the bound somewhere (else this tests nothing), this should be one: %i \n',any(crossed))
fprintf('PType constraints: (1b) and leaves it alone somewhere (else nothing checks what is NOT cut), this should be one: %i \n',any(~crossed))
fprintf('PType constraints: (1b) N_i=1 vs no-PType, one constrained update, this should be zero, PricePath r: %.3e, w: %.3e \n',max(abs(gather(PricePath2c_PT1.r-PricePath2c_NoPT.r))),max(abs(gather(PricePath2c_PT1.w-PricePath2c_NoPT.w))))
fprintf('PType constraints: (1b) N_i=2 identical types vs no-PType, one constrained update, this should be zero, PricePath r: %.3e, w: %.3e \n',max(abs(gather(PricePath2c_PT2.r-PricePath2c_NoPT.r))),max(abs(gather(PricePath2c_PT2.w-PricePath2c_NoPT.w))))
fprintf('PType constraints: (1b) N_i=2 identical types vs no-PType, one constrained update, this should be zero, GEcondnPath CapitalMarket: %.3e, LabourMarket: %.3e \n',max(abs(gather(GEcondnPath2c_PT2.CapitalMarket-GEcondnPath2c_NoPT.CapitalMarket))),max(abs(gather(GEcondnPath2c_PT2.LabourMarket-GEcondnPath2c_NoPT.LabourMarket))))
fprintf('PType constraints: (1b) N_i=2, one constrained update vs the rule applied by hand, this should be zero, r: %.3e \n',max(abs(gather(PricePath2c_PT2.r(1:T-1))-rbyhand)))
fprintf('PType constraints: (1b) N_i=2, constrained vs unconstrained update, this should NOT be zero: %.3e \n',max(abs(gather(PricePath2c_PT2.r)-gather(PricePath2u.r))))

% (1c) The same bound, 25 iterations
transpathoptions3c=transpathoptions2c;
transpathoptions3c.maxiter=25;
PricePath3c_NoPT=TransitionPath_InfHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptions3c, simoptions, vfoptions, []);
PricePath3c_PT2=TransitionPath_InfHorz_PType(PricePathFlat, ParamPathGE, T, V_final2, AgentDist_initial2, n_d,n_a_GE,n_z, 2, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params2, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions3c, simoptions, vfoptions);
fprintf('PType constraints: (1c) N_i=2 identical types vs no-PType, binding bound, 25 shooting iterations, this should be zero, PricePath r: %.3e, w: %.3e \n',max(abs(gather(PricePath3c_PT2.r-PricePath3c_NoPT.r))),max(abs(gather(PricePath3c_PT2.w-PricePath3c_NoPT.w))))
fprintf('PType constraints: (1c) N_i=2, the constrained path stays strictly inside the bound, this should be one: %i \n',all(sign(gather(PricePath3c_PT2.r(1:T-1))-bound)==sign(r0-bound)))

% (1d) A per-price step fraction (r 0.25, and w 0.9, whose bound does not bind), one update
transpathoptions2d=transpathoptions2c;
transpathoptions2d.constraintstepfraction.r=0.25;
transpathoptions2d.constraintstepfraction.w=0.9;
PricePath2d_NoPT=TransitionPath_InfHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptions2d, simoptions, vfoptions, []);
PricePath2d_PT2=TransitionPath_InfHorz_PType(PricePathFlat, ParamPathGE, T, V_final2, AgentDist_initial2, n_d,n_a_GE,n_z, 2, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params2, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions2d, simoptions, vfoptions);
rbyhand2d=gather(PricePath2u.r(1:T-1));
rbyhand2d(crossed)=r0+0.25*(bound-r0);
fprintf('PType constraints: (1d) N_i=2, per-price step fraction (r 0.25) vs no-PType, this should be zero, PricePath r: %.3e, w: %.3e \n',max(abs(gather(PricePath2d_PT2.r-PricePath2d_NoPT.r))),max(abs(gather(PricePath2d_PT2.w-PricePath2d_NoPT.w))))
fprintf('PType constraints: (1d) N_i=2, per-price step fraction vs the rule applied by hand, this should be zero, r: %.3e \n',max(abs(gather(PricePath2d_PT2.r(1:T-1))-rbyhand2d)))
fprintf('PType constraints: (1d) N_i=2, w (whose bound does not bind) is untouched, this should be zero: %.3e \n',max(abs(gather(PricePath2d_PT2.w)-gather(PricePath2u.w))))
fprintf('PType constraints: (1d) N_i=2, per-price 0.25 vs the default 0.5, this should NOT be zero: %.3e \n',max(abs(gather(PricePath2d_PT2.r-PricePath2c_PT2.r))))

%% (2) Anderson, Type-II and Type-I, 10 iterations, against a binding bound
% The bound is put halfway along the largest move of r in six unconstrained Anderson iterations without permanent types
% (as in CoreInfHorzTPath_constraints (5)). With identical types the PType Anderson takes the iterates of the one
% without permanent types, so the constrained paths must be identical.
for aa=1:2
    transpathoptionsA=transpathoptionsS;
    transpathoptionsA=rmfield(transpathoptionsA,'GEnewprice3');
    transpathoptionsA.GEnewprice=2;
    transpathoptionsA.GEnewprice2.howtoupdate=transpathoptionsS.GEnewprice3.howtoupdate;
    if aa==1
        atype='Type-II';
    else
        atype='Type-I';
        transpathoptionsA.anderson.type='I';
    end
    transpathoptionsA.maxiter=6;
    PricePath5u=TransitionPath_InfHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsA, simoptions, vfoptions, []);
    moves5=gather(PricePath5u.r(1:T-1))-r0;
    transpathoptionsAc=transpathoptionsA;
    transpathoptionsAc.maxiter=10;
    transpathoptionsAc.constrainpositive={'w'};
    transpathoptionsAc.constrainAtoB={'r'};
    if max(-moves5)>=max(moves5)
        bound5=r0-0.5*max(-moves5);
        transpathoptionsAc.constrainAtoBlimits.r=[bound5,1];
    else
        bound5=r0+0.5*max(moves5);
        transpathoptionsAc.constrainAtoBlimits.r=[-1,bound5];
    end
    transpathoptionsAu=transpathoptionsA;
    transpathoptionsAu.maxiter=10;
    PricePathAu_PT2=TransitionPath_InfHorz_PType(PricePathFlat, ParamPathGE, T, V_final2, AgentDist_initial2, n_d,n_a_GE,n_z, 2, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params2, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsAu, simoptions, vfoptions);
    [PricePathAc_NoPT,GEcondnPathAc_NoPT]=TransitionPath_InfHorz(PricePathFlat, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsAc, simoptions, vfoptions, []);
    [PricePathAc_PT1,GEcondnPathAc_PT1]=TransitionPath_InfHorz_PType(PricePathFlat, ParamPathGE, T, V_final1, AgentDist_initial1, n_d,n_a_GE,n_z, 1, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params1, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsAc, simoptions, vfoptions);
    [PricePathAc_PT2,GEcondnPathAc_PT2]=TransitionPath_InfHorz_PType(PricePathFlat, ParamPathGE, T, V_final2, AgentDist_initial2, n_d,n_a_GE,n_z, 2, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params2, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsAc, simoptions, vfoptions);
    if max(-moves5)>=max(moves5)
        crossed5=any(gather(PricePathAu_PT2.r(1:T-1))<=bound5);
    else
        crossed5=any(gather(PricePathAu_PT2.r(1:T-1))>=bound5);
    end
    fprintf('PType constraints: (2) Anderson %s, the unconstrained N_i=2 path goes on or past the bound (else this tests nothing), this should be one: %i \n',atype,crossed5)
    fprintf('PType constraints: (2) Anderson %s, N_i=1 vs no-PType, binding bound, 10 iterations, this should be zero, PricePath r: %.3e, w: %.3e \n',atype,max(abs(gather(PricePathAc_PT1.r-PricePathAc_NoPT.r))),max(abs(gather(PricePathAc_PT1.w-PricePathAc_NoPT.w))))
    fprintf('PType constraints: (2) Anderson %s, N_i=2 identical types vs no-PType, binding bound, 10 iterations, this should be zero, PricePath r: %.3e, w: %.3e \n',atype,max(abs(gather(PricePathAc_PT2.r-PricePathAc_NoPT.r))),max(abs(gather(PricePathAc_PT2.w-PricePathAc_NoPT.w))))
    fprintf('PType constraints: (2) Anderson %s, N_i=2 identical types vs no-PType, binding bound, 10 iterations, this should be zero, GEcondnPath CapitalMarket: %.3e, LabourMarket: %.3e \n',atype,max(abs(gather(GEcondnPathAc_PT2.CapitalMarket-GEcondnPathAc_NoPT.CapitalMarket))),max(abs(gather(GEcondnPathAc_PT2.LabourMarket-GEcondnPathAc_NoPT.LabourMarket))))
    fprintf('PType constraints: (2) Anderson %s, N_i=2, the constrained path ends strictly inside the bound, this should be one: %i \n',atype,all(sign(gather(PricePathAc_PT2.r(1:T-1))-bound5)==sign(r0-bound5)))
    fprintf('PType constraints: (2) Anderson %s, N_i=2, constrained vs unconstrained path, this should NOT be zero: %.3e \n',atype,max(abs(gather(PricePathAc_PT2.r-PricePathAu_PT2.r))))
end

%% (3) A price by ptype under a binding bound: two islands vs two solo constrained solves
% Each ptype has its own firm, its own r and w (prices by ptype), and its own general eqm conditions (GEptype), and
% the two differ in sigma, so the PType solve is two separate problems and must equal two solo solves exactly (for
% shooting). Both islands start from the same flat r (5% above island A's equilibrium r, so the same bound is inside
% both guesses) and their own equilibrium w. The bound on r (the same for both ptypes, as a price by ptype has) is put
% halfway along island A's largest move in one unconstrained update, so it binds for A; whether it binds for B too
% does not matter, the solo solves are the oracle for both.
Names_i={'islandA','islandB'};
sigma_i=[Params.sigma,1.5*Params.sigma];
ParamsI=Params;
ParamsI.ptweights_con=[0.5;0.5];
ParamsI.sigma=struct();
V_finalI=struct(); AgentDist_initialI=struct();
Params_solo=cell(2,1); PricePath_solo=cell(2,1);
PricePathI=struct(); PricePathI.r=struct(); PricePathI.w=struct();
for ii=1:2
    nm=Names_i{ii};
    ParamsI.sigma.(nm)=sigma_i(ii);
    Params_solo{ii}=Params;
    Params_solo{ii}.sigma=sigma_i(ii);
    p_eqm_ii=HeteroAgentStationaryEqm_InfHorz(n_d, n_a_GE, n_z, 0, pi_z, d_grid, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params_solo{ii}, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptions, vfoptions);
    Params_solo{ii}.r=p_eqm_ii.r;
    Params_solo{ii}.w=p_eqm_ii.w;
    [V_finalI.(nm),Policy_ii]=ValueFnIter_InfHorz(n_d,n_a_GE,n_z,d_grid,a_grid_GE,z_grid,pi_z,ReturnFn,Params_solo{ii},DiscountFactorParamNames,[],vfoptions);
    AgentDist_initialI.(nm)=StationaryDist_InfHorz(Policy_ii,n_d,n_a_GE,n_z,pi_z,simoptions,Params_solo{ii},[]);
    if ii==1
        rI0=1.05*p_eqm_ii.r;
    end
    PricePath_solo{ii}.r=rI0*ones(1,T);
    PricePath_solo{ii}.w=p_eqm_ii.w*ones(1,T);
    PricePathI.r.(nm)=PricePath_solo{ii}.r;
    PricePathI.w.(nm)=PricePath_solo{ii}.w;
end
AgentDist_initialI.ptweights=[0.5;0.5];
transpathoptionsI=transpathoptionsS;
transpathoptionsI.GEptype={'CapitalMarket','LabourMarket'};
transpathoptionsI.maxiter=1;
transpathoptions_solo=transpathoptionsS;
transpathoptions_solo.maxiter=1;
PricePathIu=TransitionPath_InfHorz(PricePath_solo{1}, ParamPathGE, T, V_finalI.islandA, AgentDist_initialI.islandA, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params_solo{1}, DiscountFactorParamNames, transpathoptions_solo, simoptions, vfoptions, []);
movesI=gather(PricePathIu.r(1:T-1))-rI0;
transpathoptionsIc=transpathoptionsI;
transpathoptionsIc.constrainpositive={'w'};
transpathoptionsIc.constrainAtoB={'r'};
if max(-movesI)>=max(movesI)
    boundI=rI0-0.5*max(-movesI);
    transpathoptionsIc.constrainAtoBlimits.r=[boundI,1];
    crossedI=any(gather(PricePathIu.r(1:T-1))<=boundI);
else
    boundI=rI0+0.5*max(movesI);
    transpathoptionsIc.constrainAtoBlimits.r=[-1,boundI];
    crossedI=any(gather(PricePathIu.r(1:T-1))>=boundI);
end
fprintf('PType constraints: (3) islands, the unconstrained update of island A reaches the bound (else this tests nothing), this should be one: %i \n',crossedI)
for mm=[1,10]
    transpathoptionsIc.maxiter=mm;
    transpathoptions_soloc=rmfield(transpathoptionsIc,'GEptype');
    PricePathIc=TransitionPath_InfHorz_PType(PricePathI, ParamPathGE, T, V_finalI, AgentDist_initialI, n_d,n_a_GE,n_z, Names_i, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsI, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsIc, simoptions, vfoptions);
    for ii=1:2
        nm=Names_i{ii};
        PricePath_soloc=TransitionPath_InfHorz(PricePath_solo{ii}, ParamPathGE, T, V_finalI.(nm), AgentDist_initialI.(nm), n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params_solo{ii}, DiscountFactorParamNames, transpathoptions_soloc, simoptions, vfoptions, []);
        fprintf('PType constraints: (3) islands, binding bound on r by ptype, %i shooting iteration(s), PType vs solo (island %s), this should be zero, PricePath r: %.3e, w: %.3e \n',mm,nm,max(abs(gather(PricePathIc.r.(nm)-PricePath_soloc.r))),max(abs(gather(PricePathIc.w.(nm)-PricePath_soloc.w))))
    end
    fprintf('PType constraints: (3) islands, %i shooting iteration(s), island A stays strictly inside the bound, this should be one: %i \n',mm,all(sign(gather(PricePathIc.r.islandA(1:T-1))-boundI)==sign(rI0-boundI)))
end
% Anderson with memory=0 is the shooting algorithm, so must equal it exactly, at 10 iterations
transpathoptionsIA=rmfield(transpathoptionsIc,'GEnewprice3');
transpathoptionsIA.GEnewprice=2;
transpathoptionsIA.GEnewprice2.howtoupdate=transpathoptionsIc.GEnewprice3.howtoupdate;
transpathoptionsIA.maxiter=10;
transpathoptionsIA.anderson.memory=0;
PricePathIA0=TransitionPath_InfHorz_PType(PricePathI, ParamPathGE, T, V_finalI, AgentDist_initialI, n_d,n_a_GE,n_z, Names_i, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsI, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsIA, simoptions, vfoptions);
maxdiff=0;
for ii=1:2
    nm=Names_i{ii};
    maxdiff=max([maxdiff,max(abs(gather(PricePathIA0.r.(nm)-PricePathIc.r.(nm)))),max(abs(gather(PricePathIA0.w.(nm)-PricePathIc.w.(nm))))]);
end
fprintf('PType constraints: (3) islands, binding bound, Anderson memory=0 vs shooting, 10 iterations, both islands, this should be zero: %.3e \n',maxdiff)
% And with the default memory, Anderson must stay strictly inside the bound
transpathoptionsIA=rmfield(transpathoptionsIA,'anderson');
PricePathIA5=TransitionPath_InfHorz_PType(PricePathI, ParamPathGE, T, V_finalI, AgentDist_initialI, n_d,n_a_GE,n_z, Names_i, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsI, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsIA, simoptions, vfoptions);
fprintf('PType constraints: (3) islands, binding bound, Anderson (default memory), 10 iterations, both islands strictly inside the bound, this should be one: %i \n',all(sign(gather(PricePathIA5.r.islandA(1:T-1))-boundI)==sign(rI0-boundI)) && all(sign(gather(PricePathIA5.r.islandB(1:T-1))-boundI)==sign(rI0-boundI)))
% NOT against shooting: here the bound binds in every iteration, so every Anderson step leaves the bounds, is not
% taken, and restarts the history, and Anderson then takes exactly the shooting iterates (correctly; the run of
% 2026-10-30 had shooting and Anderson bit-identical, with 4 Anderson steps refused). Against Anderson without the
% bound instead, which shows that the bound changed what Anderson did.
transpathoptionsIAu=rmfield(transpathoptionsIA,{'constrainpositive','constrainAtoB','constrainAtoBlimits'});
PricePathIA5u=TransitionPath_InfHorz_PType(PricePathI, ParamPathGE, T, V_finalI, AgentDist_initialI, n_d,n_a_GE,n_z, Names_i, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsI, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsIAu, simoptions, vfoptions);
fprintf('PType constraints: (3) islands, Anderson (default memory), constrained vs unconstrained, 10 iterations, this should NOT be zero: %.3e \n',max(abs(gather(PricePathIA5.r.islandA-PricePathIA5u.r.islandA))))

%% (4) Refusals, by name
% The eight of CoreInfHorzTPath_constraints, on the PType command (N_i=2 identical types), plus one that only a price
% by ptype can have: an initial guess that is outside its bound for one ptype only
err=zeros(1,9);
for ee=1:8
    transpathoptionsE=transpathoptionsS; transpathoptionsE.maxiter=1;
    switch ee
        case 1 % a name that is not a price in PricePath
            transpathoptionsE.constrainpositive={'sigma'}; msg='is not a price in PricePath';
        case 2 % the same price under two kinds of constraint
            transpathoptionsE.constrainpositive={'w'}; transpathoptionsE.constrain0to1={'w'}; msg='appears in more than one';
        case 3 % constrainAtoB without its limits
            transpathoptionsE.constrainAtoB={'r'}; msg='you must set transpathoptions.constrainAtoBlimits.r';
        case 4 % limits the wrong way round
            transpathoptionsE.constrainAtoB={'r'}; transpathoptionsE.constrainAtoBlimits.r=[0.2,-0.05]; msg='must be [A,B] with A<B';
        case 5 % a step fraction of one would put the price on the bound
            transpathoptionsE.constrainpositive={'w'}; transpathoptionsE.constraintstepfraction=1; msg='constraintstepfraction must be strictly between 0 and 1';
        case 6 % an initial guess outside the bound (the bumped r goes up to 5% above its equilibrium)
            transpathoptionsE.constrainAtoB={'r'}; transpathoptionsE.constrainAtoBlimits.r=[-0.05,p_eqm.r]; msg='is not strictly inside its constraint';
        case 7 % a per-price step fraction for a price that is not constrained
            transpathoptionsE.constrainpositive={'w'}; transpathoptionsE.constraintstepfraction.r=0.25; msg='is not a constrained price';
        case 8 % a step fraction that is neither a scalar nor a struct by name
            transpathoptionsE.constrainpositive={'w'}; transpathoptionsE.constraintstepfraction=[0.5,0.5]; msg='must be a scalar, or a struct with one field per constrained price';
    end
    try
        TransitionPath_InfHorz_PType(PricePathBumped, ParamPathGE, T, V_final2, AgentDist_initial2, n_d,n_a_GE,n_z, 2, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params2, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsE, simoptions, vfoptions);
    catch ME
        err(ee)=contains(ME.message,msg);
    end
end
% An initial guess outside its bound for one ptype only: island B's r guess is lowered below a lower bound that
% island A's guess is above
PricePathE=PricePathI;
PricePathE.r.islandB=0.5*rI0*ones(1,T);
transpathoptionsE=transpathoptionsI;
transpathoptionsE.constrainAtoB={'r'};
transpathoptionsE.constrainAtoBlimits.r=[0.75*rI0,1];
try
    TransitionPath_InfHorz_PType(PricePathE, ParamPathGE, T, V_finalI, AgentDist_initialI, n_d,n_a_GE,n_z, Names_i, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsI, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsE, simoptions, vfoptions);
catch ME
    err(9)=contains(ME.message,'is not strictly inside its constraint');
end
fprintf('PType constraints: (4) a constrained name that is not a price is refused by name, should be one: %i \n',err(1))
fprintf('PType constraints: (4) a price under two kinds of constraint is refused by name, should be one: %i \n',err(2))
fprintf('PType constraints: (4) constrainAtoB without its limits is refused by name, should be one: %i \n',err(3))
fprintf('PType constraints: (4) limits with A>=B are refused by name, should be one: %i \n',err(4))
fprintf('PType constraints: (4) constraintstepfraction=1 is refused by name, should be one: %i \n',err(5))
fprintf('PType constraints: (4) an initial guess outside its constraint is refused by name, should be one: %i \n',err(6))
fprintf('PType constraints: (4) a per-price step fraction for an unconstrained price is refused by name, should be one: %i \n',err(7))
fprintf('PType constraints: (4) a step fraction that is neither scalar nor struct is refused by name, should be one: %i \n',err(8))
fprintf('PType constraints: (4) an initial guess outside its constraint for one ptype only is refused by name, should be one: %i \n',err(9))

output=struct();

end
