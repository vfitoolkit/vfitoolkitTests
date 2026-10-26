function output=CoreFHorzTPathPType_Demographics_GEbequests(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,AgeWeightParamNames,transpathoptionsbaseline)
% (7) of the demographics by ptype: a general eqm path through a demographic transition that differs by
% ptype (mewj and sj both on the ParamPath, by ptype), with bequests (1-sj)*aprime received next period
% (beq_tminus1) and an interest rate, all by ptype (transpathoptions.GEptype). Each ptype is its own
% economy, so the PType solve is two solo solves without permanent types side by side. The PType
% counterpart of (9) in CoreFHorzTPath_demographics.
%   (7a) shooting at a fixed 10 iterations is separable by ptype: PType vs solo, exact
%   (7b) Anderson and shooting both reach toleranceGEcondns (flags) and agree on the eqm path
% The initial value of beq_tminus1 is by ptype (a structure). fastOLG=1.

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];
DiscountFactorParamNames={'beta','sj'};

Names_i={'typeA','typeB'};
N_i=2;
ptw=[0.6;0.4];
PTypeDistParamNames={'ptweights_demog'};

jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate.beq=@(aprime,a,z,sj) (1-sj)*aprime;
GeneralEqmEqns.rmarket=@(r,assets,beq_tminus1) r-0.05*(1+assets+beq_tminus1)^(-0.5);

%% The demographic transition of each type (as in CoreFHorzTPathPType_Demographics)
ramp=min((0:T-1)/(T/2),1);
sj_initial_A=[linspace(0.99,0.9,N_j-1),0]';
sj_final_A=[linspace(0.995,0.95,N_j-1),0]';
n_initial_A=0.01;
n_final_A=0;
sjPath_A=sj_initial_A+(sj_final_A-sj_initial_A).*ramp;
nPath_A=n_initial_A+(n_final_A-n_initial_A)*ramp;
massPath_A=zeros(N_j,T);
massPath_A(1,1)=1;
for jj=2:N_j
    massPath_A(jj,1)=massPath_A(jj-1,1)*sj_initial_A(jj-1)/(1+n_initial_A);
end
for tt=2:T
    massPath_A(1,tt)=massPath_A(1,tt-1)*(1+nPath_A(tt));
    massPath_A(2:N_j,tt)=massPath_A(1:N_j-1,tt-1).*sjPath_A(1:N_j-1,tt-1);
end
mewjPath_A=massPath_A./sum(massPath_A,1);
sj_initial_B=[linspace(0.98,0.85,N_j-1),0]';
sj_final_B=[linspace(0.99,0.92,N_j-1),0]';
n_initial_B=0.02;
n_final_B=-0.005;
sjPath_B=sj_initial_B+(sj_final_B-sj_initial_B).*ramp;
nPath_B=n_initial_B+(n_final_B-n_initial_B)*ramp;
massPath_B=zeros(N_j,T);
massPath_B(1,1)=1;
for jj=2:N_j
    massPath_B(jj,1)=massPath_B(jj-1,1)*sj_initial_B(jj-1)/(1+n_initial_B);
end
for tt=2:T
    massPath_B(1,tt)=massPath_B(1,tt-1)*(1+nPath_B(tt));
    massPath_B(2:N_j,tt)=massPath_B(1:N_j-1,tt-1).*sjPath_B(1:N_j-1,tt-1);
end
mewjPath_B=massPath_B./sum(massPath_B,1);

mewjPath_byi={mewjPath_A,mewjPath_B};
sjPath_byi={sjPath_A,sjPath_B};
sj_initial_byi={sj_initial_A,sj_initial_B};
kappa_byi={Params.kappa_j,1.3*Params.kappa_j};
beq_initial=[0.02,0.03];

%% Each type on its own
Params_solo=cell(N_i,1); V_final_PT=struct(); AgentDist_initial_PT=struct(); ParamPath_solo=cell(N_i,1);
for ii=1:N_i
    nm=Names_i{ii};
    Params_solo{ii}=Params;
    Params_solo{ii}.kappa_j=kappa_byi{ii};
    Params_solo{ii}.sj=sj_initial_byi{ii}';
    Params_solo{ii}.mewj=mewjPath_byi{ii}(:,1)';
    [V_final_PT.(nm),Policy_ii]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_solo{ii},DiscountFactorParamNames,[],vfoptions);
    AgentDist_initial_PT.(nm)=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_ii,n_d,n_a,n_z,N_j,pi_z,Params_solo{ii},simoptions);
    ParamPath_solo{ii}=ParamPath;
    ParamPath_solo{ii}.mewj=mewjPath_byi{ii};
    ParamPath_solo{ii}.sj=sjPath_byi{ii};
end
AgentDist_initial_PT.ptweights=ptw;
Params_PT=Params;
Params_PT.kappa_j=struct('typeA',kappa_byi{1},'typeB',kappa_byi{2});
Params_PT.sj=struct('typeA',sj_initial_A','typeB',sj_initial_B');
Params_PT.mewj=struct('typeA',mewjPath_A(:,1)','typeB',mewjPath_B(:,1)');
Params_PT.ptweights_demog=ptw;
ParamPath_PT=ParamPath;
ParamPath_PT.mewj=struct('typeA',mewjPath_A,'typeB',mewjPath_B);
ParamPath_PT.sj=struct('typeA',sjPath_A,'typeB',sjPath_B);
PricePath_PT.r=repmat(PricePath.r,N_i,1); % N_i-by-T, an interest rate by ptype

transpathoptionsS=transpathoptionsbaseline; % fastOLG=1
transpathoptionsS.GEnewprice=3;
transpathoptionsS.GEnewprice3.howtoupdate={'rmarket','r',0,0.5};
transpathoptionsS.verbose=0;
transpathoptionsS_PT=transpathoptionsS;
transpathoptionsS_PT.GEptype={'rmarket'};
transpathoptionsS_PT.initialvalues.beq=struct('typeA',beq_initial(1),'typeB',beq_initial(2));

%% (7a) shooting at a fixed 10 iterations, PType vs solo
transpathoptions7a=transpathoptionsS_PT;
transpathoptions7a.maxiter=10;
PricePath7a=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPath_PT, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions7a, simoptions, vfoptions);
for ii=1:N_i
    nm=Names_i{ii};
    transpathoptions7asolo=transpathoptionsS;
    transpathoptions7asolo.maxiter=10;
    transpathoptions7asolo.initialvalues.beq=beq_initial(ii);
    PricePath7asolo=TransitionPath_Case1_FHorz(PricePath, ParamPath_solo{ii}, T, V_final_PT.(nm), AgentDist_initial_PT.(nm), jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_solo{ii}, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions7asolo, simoptions, vfoptions);
    fprintf('Demographics by ptype GE path: (7a) shooting at maxiter=10, PType vs solo (type %s), this should be zero: %.3e \n',nm,max(abs(PricePath7a.r.(nm)(:)-PricePath7asolo.r(:))))
end
fprintf('Demographics by ptype GE path: (7a) for scale, 10 shooting iterations moved r by %.3e (typeA; must not be zero) \n',max(abs(PricePath7a.r.typeA(:)-PricePath.r(:))))

%% (7c) Anderson with memory=0 at the same fixed 10 iterations: every step is a plain shooting step, so == (7a)
% (Anderson is not separable by ptype, so it is not compared with solo solves; with a history it must then differ)
transpathoptions7c=rmfield(transpathoptions7a,'GEnewprice3');
transpathoptions7c.GEnewprice=2;
transpathoptions7c.GEnewprice2.howtoupdate={'rmarket','r',0,0.5};
transpathoptions7c.anderson.warmup=0;
transpathoptions7c.anderson.safeguard=0;
transpathoptions7c.anderson.memory=0;
PricePath7c0=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPath_PT, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions7c, simoptions, vfoptions);
transpathoptions7c=rmfield(transpathoptions7c,'anderson');
transpathoptions7c.anderson.warmup=0;
transpathoptions7c.anderson.safeguard=0;
PricePath7c5=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPath_PT, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions7c, simoptions, vfoptions);
fprintf('Demographics by ptype GE path: (7c) Anderson with memory=0 vs shooting, both at maxiter=10, both types, this should be zero: %.3e \n',max([max(abs(PricePath7c0.r.typeA(:)-PricePath7a.r.typeA(:))),max(abs(PricePath7c0.r.typeB(:)-PricePath7a.r.typeB(:)))]))
fprintf('Demographics by ptype GE path: (7c) Anderson with the default memory vs shooting, both at maxiter=10, both types, this should NOT be zero: %.3e \n',max([max(abs(PricePath7c5.r.typeA(:)-PricePath7a.r.typeA(:))),max(abs(PricePath7c5.r.typeB(:)-PricePath7a.r.typeB(:)))]))

%% (7b) to convergence: Anderson and shooting
transpathoptions7bS=transpathoptionsS_PT;
transpathoptions7bS.maxiter=500;
[PricePath7bS,GEcondnPath7bS]=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPath_PT, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions7bS, simoptions, vfoptions);
transpathoptions7bA=rmfield(transpathoptions7bS,'GEnewprice3');
transpathoptions7bA.GEnewprice=2;
transpathoptions7bA.GEnewprice2.howtoupdate={'rmarket','r',0,0.5};
transpathoptions7bA.maxiter=200;
[PricePath7bA,GEcondnPath7bA]=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPath_PT, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions7bA, simoptions, vfoptions);
distS=max(sqrt(GEcondnPath7bS.rmarket.typeA.^2+GEcondnPath7bS.rmarket.typeB.^2));
distA=max(sqrt(GEcondnPath7bA.rmarket.typeA.^2+GEcondnPath7bA.rmarket.typeB.^2));
fprintf('Demographics by ptype GE path: (7b) shooting reached toleranceGEcondns, should be one: %i (GE condn distance %.3e, toleranceGEcondns 1e-4) \n',distS<=1e-4,distS)
fprintf('Demographics by ptype GE path: (7b) Anderson reached toleranceGEcondns, should be one: %i (GE condn distance %.3e, toleranceGEcondns 1e-4) \n',distA<=1e-4,distA)
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('Demographics by ptype GE path: (7b) Anderson vs shooting, same eqm path (type %s), should be small: %.3e \n',nm,max(abs(PricePath7bA.r.(nm)(:)-PricePath7bS.r.(nm)(:))))
end
fprintf('Demographics by ptype GE path: (7b) for scale, the eqm r paths of the two types differ by %.3e (must not be zero) \n',max(abs(PricePath7bS.r.typeA(:)-PricePath7bS.r.typeB(:))))

output=struct();

end
