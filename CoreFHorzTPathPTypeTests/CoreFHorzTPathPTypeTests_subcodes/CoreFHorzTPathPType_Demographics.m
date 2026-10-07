function output=CoreFHorzTPathPType_Demographics(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,AgeWeightParamNames,fastOLG)
% Age weights (mewj), conditional survival (sj) and the age-1 distribution (jequaloneDist) on the
% general eqm transition path with permanent types, including paths that differ by ptype. Two types
% (kappa_j and the whole demography differ by type), on the life-cycle model of
% CoreFHorzTPath_nod_z_noe_nosemiz, with transpathoptions.fastOLG=fastOLG.
%
% THE METHOD is that of CoreFHorzTPath_demographics: the path is evaluated without moving it (shooting
% with maxiter=1 returns the conditions at the path it was handed), and the aggregates are exposed as
% the path's general eqm condition by a GE eqn that just returns them, g=assets+10*earnings. Two oracles:
%   decoupled: the GE eqn is conditional on ptype (GEptype), so it returns each ptype's own aggregates,
%              which must equal those of a solo solve without permanent types given that ptype's inputs
%   pooled:    the GE eqn is not conditional on ptype, and g is linear in the aggregates, so it must
%              equal the ptweight-weighted sum of the two solo g's
%   (1) mewj on the ParamPath, common to both types (N_j-by-T)                   decoupled
%   (2) mewj on the ParamPath, by ptype (a structure of N_j-by-T, different paths) decoupled and pooled
%   (2') as (2), with the structure's fields in the other order (typeB first): must equal (2) exactly
%   (2'') mewj on the ParamPath for typeA only (a structure with just that field): typeB is not on the
%       path for mewj, so has its constant mewj in every period                  decoupled
%   (3) mewj by ptype as a constant parameter (a structure), no path             decoupled and pooled
%   (4) (2) against constant mewj with the ratio mewj(j,t)/mewj(j,1) by ptype folded into the FnsToEvaluate
%       as an age-dependent ParamPath parameter by ptype: the same sum, so the same aggregates (pooled)
%   (5) sj on the ParamPath, by ptype, in DiscountFactorParamNames={'beta','sj'}  decoupled
%   (6) jequaloneDist as a T-path, different for each ptype                      decoupled
%   (8) AgentDistOnTransPath + EvalFnOnTransPath_AggVars (PType) with the mewj path by ptype against the
%       same standalone commands without PType, one ptype at a time
% Throughout, sj also differs by ptype as a constant (a structure), which is in the value fn.
% ((7), the general eqm path with bequests, is CoreFHorzTPathPType_Demographics_GEbequests.)

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];
DiscountFactorParamNames={'beta','sj'};

Names_i={'typeA','typeB'};
N_i=2;
ptw=[0.6;0.4];
PTypeDistParamNames={'ptweights_demog'};

% zero assets, mid point for the shock
jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;

%% The demographic transition of each type: survival rises and births change, over the first half of the path
ramp=min((0:T-1)/(T/2),1); % 0 to 1 over the first half of the path
% typeA (as in CoreFHorzTPath_demographics)
sj_initial_A=[linspace(0.99,0.9,N_j-1),0]';
sj_final_A=[linspace(0.995,0.95,N_j-1),0]';
n_initial_A=0.01;
n_final_A=0;
sjPath_A=sj_initial_A+(sj_final_A-sj_initial_A).*ramp; % N_j-by-T
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
mewjPath_A=massPath_A./sum(massPath_A,1); % N_j-by-T
% typeB: lower survival that rises more, and births that go from growing to shrinking
sj_initial_B=[linspace(0.98,0.85,N_j-1),0]';
sj_final_B=[linspace(0.99,0.92,N_j-1),0]';
n_initial_B=0.02;
n_final_B=-0.005;
sjPath_B=sj_initial_B+(sj_final_B-sj_initial_B).*ramp; % N_j-by-T
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
mewjPath_B=massPath_B./sum(massPath_B,1); % N_j-by-T
fprintf('Demographics by ptype (fastOLG=%i): the age weights move by up to %.3e (typeA) and %.3e (typeB), and differ across types by up to %.3e (none may be zero) \n',fastOLG,max(max(abs(mewjPath_A-mewjPath_A(:,1)))),max(max(abs(mewjPath_B-mewjPath_B(:,1)))),max(max(abs(mewjPath_A-mewjPath_B))))

mewjPath_byi={mewjPath_A,mewjPath_B};
sjPath_byi={sjPath_A,sjPath_B};
sj_initial_byi={sj_initial_A,sj_initial_B};
kappa_byi={Params.kappa_j,1.3*Params.kappa_j};

%% Each type on its own: Params (period-1 demographics), the value fn, V_final and Policy_final
% V_final is the stationary value fn at the period-1 demographics (as in CoreFHorzTPath_demographics)
Params_solo=cell(N_i,1); Policy_solo=cell(N_i,1); V_final_PT=struct(); Policy_final_PT=struct();
for ii=1:N_i
    Params_solo{ii}=Params;
    Params_solo{ii}.kappa_j=kappa_byi{ii};
    Params_solo{ii}.sj=sj_initial_byi{ii}';
    Params_solo{ii}.mewj=mewjPath_byi{ii}(:,1)';
    [V_final_ii,Policy_solo{ii}]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_solo{ii},DiscountFactorParamNames,[],vfoptions);
    V_final_PT.(Names_i{ii})=V_final_ii;
    Policy_final_PT.(Names_i{ii})=Policy_solo{ii};
end
% Params for the PType solves: kappa_j and sj by ptype (structures); mewj is set per check
Params_PT=Params;
Params_PT.kappa_j=struct('typeA',kappa_byi{1},'typeB',kappa_byi{2});
Params_PT.sj=struct('typeA',sj_initial_A','typeB',sj_initial_B');
Params_PT.ptweights_demog=ptw;

%% Evaluate-only: one shooting iteration returns the conditions at the path it was given
% pooled: one r, g on the pooled aggregates. decoupled: r by ptype, g conditional on ptype (GEptype).
GeneralEqmEqns.g=@(assets,earnings) assets+10*earnings;
transpathoptionsE=struct();
transpathoptionsE.fastOLG=fastOLG;
transpathoptionsE.maxiter=1;
transpathoptionsE.GEnewprice=3;
transpathoptionsE.GEnewprice3.howtoupdate={'g','r',0,0.1};
transpathoptionsE.verbose=0;
transpathoptionsED=transpathoptionsE;
transpathoptionsED.GEptype={'g'};
PricePathD.r=repmat(PricePath.r,N_i,1); % N_i-by-T, a price by ptype (needed for a GE condition by ptype)

%% (1) mewj on the ParamPath, common to both types
% Both types have the typeA age weights
ParamPath1=ParamPath;
ParamPath1.mewj=mewjPath_A;
Params1=Params_PT;
Params1.mewj=mewjPath_A(:,1)';
AgentDist_initial1=struct();
for ii=1:N_i
    Params_ii=Params_solo{ii};
    Params_ii.mewj=mewjPath_A(:,1)';
    AgentDist_initial1.(Names_i{ii})=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_solo{ii},n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions);
end
AgentDist_initial1.ptweights=ptw;
[~,GEcondnPath1]=TransitionPath_Case1_FHorz_PType(PricePathD, ParamPath1, T, V_final_PT, AgentDist_initial1, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params1, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsED, simoptions, vfoptions);
for ii=1:N_i
    Params_ii=Params_solo{ii};
    Params_ii.mewj=mewjPath_A(:,1)';
    [~,GEcondnPath1solo]=TransitionPath_Case1_FHorz(PricePath, ParamPath1, T, V_final_PT.(Names_i{ii}), AgentDist_initial1.(Names_i{ii}), jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_ii, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('Demographics by ptype (fastOLG=%i): (1) mewj path common to both types, decoupled vs solo (type %s), this should be zero: %.3e \n',fastOLG,Names_i{ii},max(abs(GEcondnPath1.g.(Names_i{ii})(:)-GEcondnPath1solo.g(:))))
end

%% (2) mewj on the ParamPath, by ptype
ParamPath2=ParamPath;
ParamPath2.mewj=struct('typeA',mewjPath_A,'typeB',mewjPath_B);
Params2=Params_PT;
Params2.mewj=struct('typeA',mewjPath_A(:,1)','typeB',mewjPath_B(:,1)');
AgentDist_initial2=struct();
g2solo=cell(N_i,1);
for ii=1:N_i
    AgentDist_initial2.(Names_i{ii})=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_solo{ii},n_d,n_a,n_z,N_j,pi_z,Params_solo{ii},simoptions);
end
AgentDist_initial2.ptweights=ptw;
[~,GEcondnPath2D]=TransitionPath_Case1_FHorz_PType(PricePathD, ParamPath2, T, V_final_PT, AgentDist_initial2, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params2, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsED, simoptions, vfoptions);
[~,GEcondnPath2P]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath2, T, V_final_PT, AgentDist_initial2, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params2, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsE, simoptions, vfoptions);
for ii=1:N_i
    ParamPath_ii=ParamPath;
    ParamPath_ii.mewj=mewjPath_byi{ii};
    [~,GEcondnPath2solo]=TransitionPath_Case1_FHorz(PricePath, ParamPath_ii, T, V_final_PT.(Names_i{ii}), AgentDist_initial2.(Names_i{ii}), jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_solo{ii}, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    g2solo{ii}=GEcondnPath2solo.g;
    fprintf('Demographics by ptype (fastOLG=%i): (2) mewj path by ptype, decoupled vs solo (type %s), this should be zero: %.3e \n',fastOLG,Names_i{ii},max(abs(GEcondnPath2D.g.(Names_i{ii})(:)-g2solo{ii}(:))))
end
fprintf('Demographics by ptype (fastOLG=%i): (2) mewj path by ptype, pooled vs ptweight-weighted solo, this should be zero: %.3e \n',fastOLG,max(abs(GEcondnPath2P.g(:)-(ptw(1)*g2solo{1}(:)+ptw(2)*g2solo{2}(:)))))

%% (2') the same, with the fields of the ParamPath structure in the other order
ParamPath2r=ParamPath;
ParamPath2r.mewj=struct('typeB',mewjPath_B,'typeA',mewjPath_A);
[~,GEcondnPath2rD]=TransitionPath_Case1_FHorz_PType(PricePathD, ParamPath2r, T, V_final_PT, AgentDist_initial2, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params2, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsED, simoptions, vfoptions);
for ii=1:N_i
    fprintf('Demographics by ptype (fastOLG=%i): (2'') mewj path by ptype with the fields typeB first vs (2) (type %s), this should be zero: %.3e \n',fastOLG,Names_i{ii},max(abs(GEcondnPath2rD.g.(Names_i{ii})(:)-GEcondnPath2D.g.(Names_i{ii})(:))))
end

%% (2'') mewj on the ParamPath for typeA only; typeB keeps its constant (period-1) mewj
ParamPath2a=ParamPath;
ParamPath2a.mewj=struct('typeA',mewjPath_A);
[~,GEcondnPath2aD]=TransitionPath_Case1_FHorz_PType(PricePathD, ParamPath2a, T, V_final_PT, AgentDist_initial2, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params2, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsED, simoptions, vfoptions);
fprintf('Demographics by ptype (fastOLG=%i): (2'''') mewj path for typeA only, decoupled vs solo with the path (type typeA), this should be zero: %.3e \n',fastOLG,max(abs(GEcondnPath2aD.g.typeA(:)-g2solo{1}(:))))
[~,GEcondnPath2asolo]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final_PT.typeB, AgentDist_initial2.typeB, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_solo{2}, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
fprintf('Demographics by ptype (fastOLG=%i): (2'''') mewj path for typeA only, decoupled vs solo without a mewj path (type typeB), this should be zero: %.3e \n',fastOLG,max(abs(GEcondnPath2aD.g.typeB(:)-GEcondnPath2asolo.g(:))))

%% (3) mewj by ptype as a constant parameter, no path
Params3=Params_PT;
Params3.mewj=struct('typeA',mewjPath_A(:,end)','typeB',mewjPath_B(:,end)'); % the final age weights, which are far apart
AgentDist_initial3=struct();
g3solo=cell(N_i,1);
for ii=1:N_i
    Params_ii=Params_solo{ii};
    Params_ii.mewj=mewjPath_byi{ii}(:,end)';
    AgentDist_initial3.(Names_i{ii})=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_solo{ii},n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions);
    [~,GEcondnPath3solo]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final_PT.(Names_i{ii}), AgentDist_initial3.(Names_i{ii}), jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_ii, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    g3solo{ii}=GEcondnPath3solo.g;
end
AgentDist_initial3.ptweights=ptw;
[~,GEcondnPath3D]=TransitionPath_Case1_FHorz_PType(PricePathD, ParamPath, T, V_final_PT, AgentDist_initial3, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params3, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsED, simoptions, vfoptions);
[~,GEcondnPath3P]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial3, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params3, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsE, simoptions, vfoptions);
for ii=1:N_i
    fprintf('Demographics by ptype (fastOLG=%i): (3) constant mewj by ptype, decoupled vs solo (type %s), this should be zero: %.3e \n',fastOLG,Names_i{ii},max(abs(GEcondnPath3D.g.(Names_i{ii})(:)-g3solo{ii}(:))))
end
fprintf('Demographics by ptype (fastOLG=%i): (3) constant mewj by ptype, pooled vs ptweight-weighted solo, this should be zero: %.3e \n',fastOLG,max(abs(GEcondnPath3P.g(:)-(ptw(1)*g3solo{1}(:)+ptw(2)*g3solo{2}(:)))))

%% (4) (2) against constant mewj with the ratio folded into the FnsToEvaluate, by ptype
% sum_j mewj(j,t)*E[f|j,t] = sum_j mewj(j,1)*E[f*mewj(j,t)/mewj(j,1)|j,t], so the aggregates agree
FnsToEvaluate4.assets=@(aprime,a,z,mewjratio) a*mewjratio;
FnsToEvaluate4.earnings=@(aprime,a,z,w,kappa_j,mewjratio) w*kappa_j*z*mewjratio;
ParamPath4=ParamPath;
ParamPath4.mewjratio=struct('typeA',mewjPath_A./mewjPath_A(:,1),'typeB',mewjPath_B./mewjPath_B(:,1)); % N_j-by-T by ptype
Params4=Params2; % period-1 mewj by ptype, held constant
Params4.mewjratio=struct('typeA',ones(1,N_j),'typeB',ones(1,N_j));
[~,GEcondnPath4P]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath4, T, V_final_PT, AgentDist_initial2, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate4, GeneralEqmEqns, Params4, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsE, simoptions, vfoptions);
fprintf('Demographics by ptype (fastOLG=%i): (4) mewj path by ptype vs constant mewj with the ratio by ptype in the FnsToEvaluate, pooled, this should be zero: %.3e \n',fastOLG,max(abs(GEcondnPath2P.g(:)-GEcondnPath4P.g(:))))

%% (5) sj on the ParamPath, by ptype (and mewj on the ParamPath by ptype, as in (2))
ParamPath5=ParamPath2;
ParamPath5.sj=struct('typeA',sjPath_A,'typeB',sjPath_B);
[~,GEcondnPath5D]=TransitionPath_Case1_FHorz_PType(PricePathD, ParamPath5, T, V_final_PT, AgentDist_initial2, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params2, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsED, simoptions, vfoptions);
for ii=1:N_i
    ParamPath_ii=ParamPath;
    ParamPath_ii.mewj=mewjPath_byi{ii};
    ParamPath_ii.sj=sjPath_byi{ii};
    [~,GEcondnPath5solo]=TransitionPath_Case1_FHorz(PricePath, ParamPath_ii, T, V_final_PT.(Names_i{ii}), AgentDist_initial2.(Names_i{ii}), jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_solo{ii}, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('Demographics by ptype (fastOLG=%i): (5) sj and mewj paths by ptype, decoupled vs solo (type %s), this should be zero: %.3e \n',fastOLG,Names_i{ii},max(abs(GEcondnPath5D.g.(Names_i{ii})(:)-GEcondnPath5solo.g(:))))
end
fprintf('Demographics by ptype (fastOLG=%i): (5) for scale, the sj path changes g by %.3e (typeA; must not be zero, or (5) adds nothing to (2)) \n',fastOLG,max(abs(GEcondnPath5D.g.typeA(:)-GEcondnPath2D.g.typeA(:))))

%% (6) jequaloneDist as a T-path, different for each ptype
% typeA newborns get more assets from mid-path, typeB from a third of the way through, at another grid point
jequaloneDist6=struct();
jequaloneDist6.typeA=zeros(n_a,n_z,T,'gpuArray');
jequaloneDist6.typeA(1,ceil(n_z/2),1:floor(T/2))=1;
jequaloneDist6.typeA(5,ceil(n_z/2),floor(T/2)+1:T)=1;
jequaloneDist6.typeB=zeros(n_a,n_z,T,'gpuArray');
jequaloneDist6.typeB(1,ceil(n_z/2),1:floor(T/3))=1;
jequaloneDist6.typeB(9,1,floor(T/3)+1:T)=1;
[~,GEcondnPath6D]=TransitionPath_Case1_FHorz_PType(PricePathD, ParamPath2, T, V_final_PT, AgentDist_initial2, jequaloneDist6, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params2, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsED, simoptions, vfoptions);
for ii=1:N_i
    ParamPath_ii=ParamPath;
    ParamPath_ii.mewj=mewjPath_byi{ii};
    [~,GEcondnPath6solo]=TransitionPath_Case1_FHorz(PricePath, ParamPath_ii, T, V_final_PT.(Names_i{ii}), AgentDist_initial2.(Names_i{ii}), jequaloneDist6.(Names_i{ii}), n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_solo{ii}, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('Demographics by ptype (fastOLG=%i): (6) jequaloneDist path by ptype, decoupled vs solo (type %s), this should be zero: %.3e \n',fastOLG,Names_i{ii},max(abs(GEcondnPath6D.g.(Names_i{ii})(:)-GEcondnPath6solo.g(:))))
end
fprintf('Demographics by ptype (fastOLG=%i): (6) for scale, the newborn paths change g by %.3e (typeB; must not be zero, or (6) adds nothing to (2)) \n',fastOLG,max(abs(GEcondnPath6D.g.typeB(:)-GEcondnPath2D.g.typeB(:))))

%% (8) The standalone commands with the mewj path by ptype, against the same commands without PType
transpathoptions8=struct();
transpathoptions8.fastOLG=fastOLG;
[~,PolicyPath8]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath2, T, V_final_PT, Policy_final_PT, Params2, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions8, vfoptions);
AgentDistPath8=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial2, jequaloneDist, PricePath, ParamPath2, PolicyPath8, AgeWeightParamNames, PTypeDistParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z, T,Params2, transpathoptions8, simoptions);
AggVarsPath8=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath8, PolicyPath8, PricePath, ParamPath2, Params2, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptions8, simoptions);
for ii=1:N_i
    nm=Names_i{ii};
    ParamPath_ii=ParamPath;
    ParamPath_ii.mewj=mewjPath_byi{ii};
    [~,PolicyPath8solo]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath_ii, T, V_final_PT.(nm), Policy_final_PT.(nm), Params_solo{ii}, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions8, vfoptions);
    AgentDistPath8solo=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial2.(nm), jequaloneDist, PricePath, ParamPath_ii, PolicyPath8solo, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params_solo{ii}, transpathoptions8, simoptions);
    AggVarsPath8solo=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath8solo, PolicyPath8solo, PricePath, ParamPath_ii, Params_solo{ii}, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptions8, simoptions);
    fprintf('Demographics by ptype (fastOLG=%i): (8) standalone commands, mewj path by ptype, PType vs solo, assets (type %s), this should be zero: %.3e \n',fastOLG,nm,max(abs(AggVarsPath8.assets.(nm).Mean(:)-AggVarsPath8solo.assets.Mean(:))))
    fprintf('Demographics by ptype (fastOLG=%i): (8) standalone commands, mewj path by ptype, PType vs solo, earnings (type %s), this should be zero: %.3e \n',fastOLG,nm,max(abs(AggVarsPath8.earnings.(nm).Mean(:)-AggVarsPath8solo.earnings.Mean(:))))
    % and the standalone aggregates against the ones inside the path (the decoupled GE condition of (2))
    fprintf('Demographics by ptype (fastOLG=%i): (8) standalone vs inside the path, mewj path by ptype (type %s), this should be zero: %.3e \n',fastOLG,nm,max(abs(gather(AggVarsPath8.assets.(nm).Mean(1:T-1)+10*AggVarsPath8.earnings.(nm).Mean(1:T-1))-gather(GEcondnPath2D.g.(nm)))))
end

output=struct();

end
