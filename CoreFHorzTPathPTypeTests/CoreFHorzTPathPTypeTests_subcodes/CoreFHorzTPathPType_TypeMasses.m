function output=CoreFHorzTPathPType_TypeMasses(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,AgeWeightParamNames,transpathoptionsbaseline,fastOLG)
% The mass of each permanent type changing over the transition path: the ptype weights
% (PTypeDistParamNames) on the ParamPath, as N_i-by-T, T-by-N_i, or a structure by ptype name. Two types
% (kappa_j differs, given as a structure), on the life-cycle model of CoreFHorzTPath_nod_z_noe_nosemiz,
% with transpathoptions.fastOLG=fastOLG.
%
% Each ptype's agent distribution has mass one, so the ptype weights only enter when the aggregates are
% pooled. THE METHOD is that of CoreFHorzTPathPType_Demographics: the path is evaluated without moving
% it (shooting with maxiter=1), and the aggregates are exposed by the GE eqn g=assets+10*earnings, which
% is linear in them. So pooled g must equal sum_i ptweights(i,t)*g_i(t), with g_i from a solo solve.
%   (1) the weights on the ParamPath but constant == the weights not on the ParamPath
%   (2) a weights path: pooled == ptweight-weighted solo g, period by period; given as N_i-by-T, as
%       T-by-N_i, and as a structure (with its fields in the opposite order to Names_i), all the same
%   (3) GE conditions by ptype (GEptype) are unaffected by the weights path: == solo
%   (4) the standalone commands (AgentDistOnTransPath, EvalFnOnTransPath_AggVars PType) with the weights
%       path == inside the path; and without it, the constant weights of (1)
%   (5) the weights path together with mewj paths by ptype: pooled == weighted solo with their own mewj
%   (6) a general eqm path with one interest rate and the weights path: shooting and Anderson converge
%       and agree (fastOLG=1 only)
%   (7) Parameters weights that differ from AgentDist_initial.ptweights (no path): the path uses
%       AgentDist_initial.ptweights, so == (1); and a path whose period 1 differs from them: period 1 uses
%       the path. Both also print a warning, which is not a check.

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];
DiscountFactorParamNames={'beta'};

Names_i={'typeA','typeB'};
N_i=2;
PTypeDistParamNames={'ptweights_mass'};

jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;

%% The type masses: 0.6/0.4 at the start, 0.35/0.65 from the middle of the path on
ptw_initial=[0.6;0.4];
ramp=min((0:T-1)/(T/2),1); % 0 to 1 over the first half of the path
ptwPath=[0.6+(0.35-0.6)*ramp; 0.4+(0.65-0.4)*ramp]; % N_i-by-T, columns sum to one

%% Each type on its own: V_final and the initial agent distribution, at the period-1 prices
kappa_byi={Params.kappa_j,1.3*Params.kappa_j};
Params_solo=cell(N_i,1); V_final_PT=struct(); Policy_final_PT=struct(); AgentDist_initial_PT=struct(); Policy_solo=cell(N_i,1);
for ii=1:N_i
    nm=Names_i{ii};
    Params_solo{ii}=Params;
    Params_solo{ii}.kappa_j=kappa_byi{ii};
    [V_final_PT.(nm),Policy_solo{ii}]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_solo{ii},DiscountFactorParamNames,[],vfoptions);
    Policy_final_PT.(nm)=Policy_solo{ii};
    AgentDist_initial_PT.(nm)=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_solo{ii},n_d,n_a,n_z,N_j,pi_z,Params_solo{ii},simoptions);
end
AgentDist_initial_PT.ptweights=ptw_initial;
Params_PT=Params;
Params_PT.kappa_j=struct('typeA',kappa_byi{1},'typeB',kappa_byi{2});
Params_PT.ptweights_mass=ptw_initial;

%% Evaluate-only: one shooting iteration returns the conditions at the path it was given
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

% The solo g of each type
gsolo=zeros(N_i,T-1);
for ii=1:N_i
    [~,GEcondnPath_solo]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final_PT.(Names_i{ii}), AgentDist_initial_PT.(Names_i{ii}), jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_solo{ii}, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    gsolo(ii,:)=gather(GEcondnPath_solo.g(:))';
end

%% (1) the weights on the ParamPath but constant == not on the ParamPath
[~,GEcondnPath1]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsE, simoptions, vfoptions);
ParamPath1c=ParamPath;
ParamPath1c.ptweights_mass=ptw_initial.*ones(1,T);
[~,GEcondnPath1c]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath1c, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsE, simoptions, vfoptions);
fprintf('Type masses (fastOLG=%i): (1) weights on the ParamPath but constant vs not on it, pooled, this should be zero: %.3e \n',fastOLG,max(abs(GEcondnPath1c.g(:)-GEcondnPath1.g(:))))
fprintf('Type masses (fastOLG=%i): (1) weights not on the ParamPath, pooled vs ptweight-weighted solo, this should be zero: %.3e \n',fastOLG,max(abs(gather(GEcondnPath1.g(:))'-sum(ptw_initial.*gsolo,1))))

%% (2) a weights path, given three ways
ParamPath2=ParamPath;
ParamPath2.ptweights_mass=ptwPath; % N_i-by-T
[~,GEcondnPath2]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath2, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsE, simoptions, vfoptions);
fprintf('Type masses (fastOLG=%i): (2) weights path N_i-by-T, pooled vs ptweight-weighted solo period by period, this should be zero: %.3e \n',fastOLG,max(abs(gather(GEcondnPath2.g(:))'-sum(ptwPath(:,1:T-1).*gsolo,1))))
fprintf('Type masses (fastOLG=%i): (2) for scale, the weights path changes pooled g by %.3e (must not be zero) \n',fastOLG,max(abs(GEcondnPath2.g(:)-GEcondnPath1.g(:))))
ParamPath2t=ParamPath;
ParamPath2t.ptweights_mass=ptwPath'; % T-by-N_i
[~,GEcondnPath2t]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath2t, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsE, simoptions, vfoptions);
fprintf('Type masses (fastOLG=%i): (2) weights path T-by-N_i vs N_i-by-T, this should be zero: %.3e \n',fastOLG,max(abs(GEcondnPath2t.g(:)-GEcondnPath2.g(:))))
ParamPath2s=ParamPath;
ParamPath2s.ptweights_mass=struct('typeB',ptwPath(2,:),'typeA',ptwPath(1,:)); % fields in the opposite order to Names_i
[~,GEcondnPath2s]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath2s, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsE, simoptions, vfoptions);
fprintf('Type masses (fastOLG=%i): (2) weights path as a structure (typeB first) vs N_i-by-T, this should be zero: %.3e \n',fastOLG,max(abs(GEcondnPath2s.g(:)-GEcondnPath2.g(:))))

%% (3) GE conditions by ptype are unaffected by the weights path
[~,GEcondnPath3]=TransitionPath_Case1_FHorz_PType(PricePathD, ParamPath2, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsED, simoptions, vfoptions);
for ii=1:N_i
    fprintf('Type masses (fastOLG=%i): (3) weights path, GE condition by ptype vs solo (type %s), this should be zero: %.3e \n',fastOLG,Names_i{ii},max(abs(gather(GEcondnPath3.g.(Names_i{ii})(:))'-gsolo(ii,:))))
end

%% (4) the standalone commands, with and without the weights path
transpathoptions4=struct();
transpathoptions4.fastOLG=fastOLG;
[~,PolicyPath4]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath2, T, V_final_PT, Policy_final_PT, Params_PT, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions4, vfoptions);
AgentDistPath4=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist, PricePath, ParamPath2, PolicyPath4, AgeWeightParamNames, PTypeDistParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z, T,Params_PT, transpathoptions4, simoptions);
AggVarsPath4=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath4, PolicyPath4, PricePath, ParamPath2, Params_PT, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptions4, simoptions);
fprintf('Type masses (fastOLG=%i): (4) AgentDistPath.ptweights with the weights path is N_i-by-T and equals it, this should be zero: %.3e \n',fastOLG,max(abs(reshape(AgentDistPath4.ptweights,[],1)-reshape(ptwPath,[],1)))+abs(numel(AgentDistPath4.ptweights)-N_i*T))
fprintf('Type masses (fastOLG=%i): (4) standalone vs inside the path, weights path, pooled, this should be zero: %.3e \n',fastOLG,max(abs(gather(AggVarsPath4.assets.Mean(1:T-1)+10*AggVarsPath4.earnings.Mean(1:T-1))-gather(GEcondnPath2.g))))
[~,PolicyPath4n]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, Policy_final_PT, Params_PT, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions4, vfoptions);
AgentDistPath4n=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist, PricePath, ParamPath, PolicyPath4n, AgeWeightParamNames, PTypeDistParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z, T,Params_PT, transpathoptions4, simoptions);
AggVarsPath4n=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath4n, PolicyPath4n, PricePath, ParamPath, Params_PT, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptions4, simoptions);
fprintf('Type masses (fastOLG=%i): (4) without the weights path, AgentDistPath.ptweights is AgentDist_initial.ptweights, this should be zero: %.3e \n',fastOLG,max(abs(reshape(AgentDistPath4n.ptweights,[],1)-ptw_initial))+abs(numel(AgentDistPath4n.ptweights)-N_i))
fprintf('Type masses (fastOLG=%i): (4) standalone vs inside the path, no weights path, pooled, this should be zero: %.3e \n',fastOLG,max(abs(gather(AggVarsPath4n.assets.Mean(1:T-1)+10*AggVarsPath4n.earnings.Mean(1:T-1))-gather(GEcondnPath1.g))))

%% (5) the weights path together with mewj paths by ptype
% typeA's age weights tilt towards the old, typeB's towards the young, both starting from Params.mewj
tilt=linspace(-1,1,N_j)';
mewjPath_A=Params.mewj(:).*(1+0.3*tilt*ramp); mewjPath_A=mewjPath_A./sum(mewjPath_A,1);
mewjPath_B=Params.mewj(:).*(1-0.3*tilt*ramp); mewjPath_B=mewjPath_B./sum(mewjPath_B,1);
ParamPath5=ParamPath2;
ParamPath5.mewj=struct('typeA',mewjPath_A,'typeB',mewjPath_B);
[~,GEcondnPath5]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath5, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsE, simoptions, vfoptions);
gsolo5=zeros(N_i,T-1);
mewjPath_byi={mewjPath_A,mewjPath_B};
for ii=1:N_i
    ParamPath_ii=ParamPath;
    ParamPath_ii.mewj=mewjPath_byi{ii};
    [~,GEcondnPath5solo]=TransitionPath_Case1_FHorz(PricePath, ParamPath_ii, T, V_final_PT.(Names_i{ii}), AgentDist_initial_PT.(Names_i{ii}), jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_solo{ii}, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    gsolo5(ii,:)=gather(GEcondnPath5solo.g(:))';
end
fprintf('Type masses (fastOLG=%i): (5) weights path and mewj paths by ptype, pooled vs ptweight-weighted solo, this should be zero: %.3e \n',fastOLG,max(abs(gather(GEcondnPath5.g(:))'-sum(ptwPath(:,1:T-1).*gsolo5,1))))
fprintf('Type masses (fastOLG=%i): (5) for scale, the mewj paths change pooled g by %.3e (must not be zero) \n',fastOLG,max(abs(GEcondnPath5.g(:)-GEcondnPath2.g(:))))

%% (6) a general eqm path with one interest rate and the weights path
if fastOLG==1
    GeneralEqmEqns6.rmarket=@(r,assets) r-0.05*(1+assets)^(-0.5);
    transpathoptions6S=transpathoptionsbaseline; % fastOLG=1
    transpathoptions6S.GEnewprice=3;
    transpathoptions6S.GEnewprice3.howtoupdate={'rmarket','r',0,0.5};
    transpathoptions6S.maxiter=500;
    transpathoptions6S.verbose=0;
    [PricePath6S,GEcondnPath6S]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath2, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns6, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions6S, simoptions, vfoptions);
    transpathoptions6A=rmfield(transpathoptions6S,'GEnewprice3');
    transpathoptions6A.GEnewprice=2;
    transpathoptions6A.GEnewprice2.howtoupdate={'rmarket','r',0,0.5};
    transpathoptions6A.maxiter=200;
    [PricePath6A,GEcondnPath6A]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath2, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns6, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions6A, simoptions, vfoptions);
    distS=max(abs(GEcondnPath6S.rmarket));
    distA=max(abs(GEcondnPath6A.rmarket));
    fprintf('Type masses GE path: (6) shooting reached toleranceGEcondns, should be one: %i (GE condn distance %.3e, toleranceGEcondns 1e-4) \n',distS<=1e-4,distS)
    fprintf('Type masses GE path: (6) Anderson reached toleranceGEcondns, should be one: %i (GE condn distance %.3e, toleranceGEcondns 1e-4) \n',distA<=1e-4,distA)
    fprintf('Type masses GE path: (6) Anderson vs shooting, same eqm path, should be small: %.3e \n',max(abs(PricePath6A.r(:)-PricePath6S.r(:))))
    fprintf('Type masses GE path: (6) for scale, the eqm r moves along the path by %.3e (must not be zero) \n',max(PricePath6S.r(:))-min(PricePath6S.r(:)))
end

%% (7) weights that differ from AgentDist_initial.ptweights (each also prints a warning, which is not a check)
% Not on the path: the path uses AgentDist_initial.ptweights, whatever Parameters says
Params7=Params_PT;
Params7.ptweights_mass=[0.5;0.5];
[~,GEcondnPath7]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params7, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsE, simoptions, vfoptions);
fprintf('Type masses (fastOLG=%i): (7) Parameters weights differ from AgentDist_initial.ptweights, the path uses the latter, this should be zero: %.3e \n',fastOLG,max(abs(GEcondnPath7.g(:)-GEcondnPath1.g(:))))
% On the path with a different period 1: period 1 uses the path
ptwPath7=ptwPath;
ptwPath7(:,1)=[0.5;0.5];
ParamPath7=ParamPath;
ParamPath7.ptweights_mass=ptwPath7;
[~,GEcondnPath7p]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath7, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsE, simoptions, vfoptions);
fprintf('Type masses (fastOLG=%i): (7) weights path with period 1 different from AgentDist_initial.ptweights, pooled vs path-weighted solo, this should be zero: %.3e \n',fastOLG,max(abs(gather(GEcondnPath7p.g(:))'-sum(ptwPath7(:,1:T-1).*gsolo,1))))

output=struct();

end
