function output=CoreFHorzTPathPType_GE_DecoupledTypes(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,fastOLG)
% The general eqm transition path with permanent types, for four types that differ in nearly
% everything a type can differ in, and whose general eqm condition is conditional on ptype
% (transpathoptions.GEptype), with an interest rate per ptype. The types then do not interact at all,
% so the one PType solve must give exactly what four separate solves without permanent types give.
%   none: no z, no e; its own n_a and a_grid; its own age weights (mewj) and survival (sj)
%   z:    z; divide-and-conquer; the age-1 distribution is a path over t
%   e:    e; a decision variable d
%   ze:   z and e, with a z_grid different from type z (same n_z, other values)
% Each solve runs the shooting algorithm for maxiter=2 (so returns the price path after one update,
% and the general eqm conditions on it), with transpathoptions.fastOLG=fastOLG. Checked: PricePath
% and GEcondnPath of each type against its solo solve, and PricePathptype_vectoroutput=1 against =0.
% The exogenous shock setup writes each type's z and e paths into transpathoptions, so types with
% different shocks are what would show a shared transpathoptions.

Names_i={'none','z','e','ze'};
N_i=length(Names_i);

n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;
n_a_none=61;
a_grid_none=5*linspace(0,1,n_a_none)'.^3;

%% Per-type inputs, as structures
n_d_PT=struct(); d_grid_PT=struct(); n_a_PT=struct(); a_grid_PT=struct();
n_z_PT=struct(); z_grid_PT=struct(); pi_z_PT=struct();
ReturnFn_PT=struct(); jequaloneDist_PT=struct();
vfopt_PT=struct(); simopt_PT=struct();
vfopt_PT.n_e=struct(); vfopt_PT.e_grid=struct(); vfopt_PT.pi_e=struct();
simopt_PT.n_e=struct(); simopt_PT.e_grid=struct(); simopt_PT.pi_e=struct();
vfopt_PT.divideandconquer=struct();
Params_PT=Params;
Params_PT.mewj=struct();
Params_PT.sj=struct();
DiscountFactorParamNames_GE={DiscountFactorParamNames{:},'sj'};

% type 1: none
nm='none';
n_d_PT.(nm)=0; d_grid_PT.(nm)=[];
n_a_PT.(nm)=n_a_none; a_grid_PT.(nm)=a_grid_none;
n_z_PT.(nm)=0; z_grid_PT.(nm)=[]; pi_z_PT.(nm)=[];
ReturnFn_PT.(nm)=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate.assets.(nm)=@(aprime,a) a;
FnsToEvaluate.earnings.(nm)=@(aprime,a,w,kappa_j) w*kappa_j;
jequaloneDist_PT.(nm)=zeros(n_a_none,1,'gpuArray');
jequaloneDist_PT.(nm)(1)=1;
Params_PT.mewj.(nm)=(0.97.^(0:1:N_j-1))/sum(0.97.^(0:1:N_j-1));
Params_PT.sj.(nm)=[0.98*ones(1,N_j-1),0];

% type 2: z
nm='z';
n_d_PT.(nm)=0; d_grid_PT.(nm)=[];
n_a_PT.(nm)=n_a; a_grid_PT.(nm)=a_grid;
n_z_PT.(nm)=n_z; z_grid_PT.(nm)=z_grid; pi_z_PT.(nm)=pi_z;
ReturnFn_PT.(nm)=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate.assets.(nm)=@(aprime,a,z) a;
FnsToEvaluate.earnings.(nm)=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
% age-1 distribution is a path: the newborns start with more assets in the second half of the path
jequaloneDist_PT.(nm)=zeros(n_a,n_z,T,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_z/2),1:floor(T/2))=1;
jequaloneDist_PT.(nm)(5,ceil(n_z/2),floor(T/2)+1:T)=1;
vfopt_PT.divideandconquer.(nm)=1;
Params_PT.mewj.(nm)=Params.mewj;
Params_PT.sj.(nm)=ones(1,N_j);

% type 3: e, with d
nm='e';
n_d_PT.(nm)=n_d; d_grid_PT.(nm)=d_grid;
n_a_PT.(nm)=n_a; a_grid_PT.(nm)=a_grid;
n_z_PT.(nm)=0; z_grid_PT.(nm)=[]; pi_z_PT.(nm)=[];
ReturnFn_PT.(nm)=@(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
FnsToEvaluate.assets.(nm)=@(d,aprime,a,e) a;
FnsToEvaluate.earnings.(nm)=@(d,aprime,a,e,w,kappa_j) w*kappa_j*d*e;
vfopt_PT.n_e.(nm)=n_e; vfopt_PT.e_grid.(nm)=e_grid; vfopt_PT.pi_e.(nm)=pi_e;
simopt_PT.n_e.(nm)=n_e; simopt_PT.e_grid.(nm)=e_grid; simopt_PT.pi_e.(nm)=pi_e;
jequaloneDist_PT.(nm)=zeros(n_a,n_e,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_e/2))=1;
Params_PT.mewj.(nm)=Params.mewj;
Params_PT.sj.(nm)=ones(1,N_j);

% type 4: z and e, z_grid differs from type z
nm='ze';
n_d_PT.(nm)=0; d_grid_PT.(nm)=[];
n_a_PT.(nm)=n_a; a_grid_PT.(nm)=a_grid;
n_z_PT.(nm)=n_z; z_grid_PT.(nm)=z_grid.^1.5; pi_z_PT.(nm)=pi_z;
ReturnFn_PT.(nm)=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate.assets.(nm)=@(aprime,a,z,e) a;
FnsToEvaluate.earnings.(nm)=@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;
vfopt_PT.n_e.(nm)=n_e; vfopt_PT.e_grid.(nm)=e_grid; vfopt_PT.pi_e.(nm)=pi_e;
simopt_PT.n_e.(nm)=n_e; simopt_PT.e_grid.(nm)=e_grid; simopt_PT.pi_e.(nm)=pi_e;
jequaloneDist_PT.(nm)=zeros(n_a,n_z,n_e,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_z/2),ceil(n_e/2))=1;
Params_PT.mewj.(nm)=Params.mewj;
Params_PT.sj.(nm)=ones(1,N_j);

%% General eqm: one condition per ptype, each with its own interest rate
GeneralEqmEqns.CapitalMarket=@(r,assets,earnings) r-0.02-0.01*assets/(1+earnings);
PricePath_PT=struct();
PricePath_PT.r=repmat(PricePath.r,N_i,1); % N_i-by-T, a price that depends on ptype

transpathoptions=struct();
transpathoptions.fastOLG=fastOLG;
transpathoptions.GEnewprice=3; % the shooting algorithm
transpathoptions.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5};
transpathoptions.GEptype={'CapitalMarket'};
transpathoptions.maxiter=2; % returns the price path after one update, and the general eqm conditions evaluated on it
transpathoptions.verbose=0;

%% Each type: the initial agent distribution (stationary at the period-1 prices) and V_final (stationary at the final prices)
V_final_PT=struct(); AgentDist_initial_PT=struct();
vfopt_solo=cell(N_i,1); simopt_solo=cell(N_i,1); Params_solo=cell(N_i,1);
for ii=1:N_i
    nm=Names_i{ii};
    vfopt_solo{ii}=struct();
    simopt_solo{ii}=struct();
    if isfield(vfopt_PT.n_e,nm)
        vfopt_solo{ii}.n_e=vfopt_PT.n_e.(nm);
        vfopt_solo{ii}.e_grid=vfopt_PT.e_grid.(nm);
        vfopt_solo{ii}.pi_e=vfopt_PT.pi_e.(nm);
        simopt_solo{ii}.n_e=vfopt_solo{ii}.n_e;
        simopt_solo{ii}.e_grid=vfopt_solo{ii}.e_grid;
        simopt_solo{ii}.pi_e=vfopt_solo{ii}.pi_e;
    end
    if isfield(vfopt_PT.divideandconquer,nm)
        vfopt_solo{ii}.divideandconquer=vfopt_PT.divideandconquer.(nm);
    end
    Params_solo{ii}=Params;
    Params_solo{ii}.mewj=Params_PT.mewj.(nm);
    Params_solo{ii}.sj=Params_PT.sj.(nm);

    % the stationary age-1 distribution is that of period 1
    if strcmp(nm,'z')
        jequaloneDist_1=jequaloneDist_PT.(nm)(:,:,1);
    else
        jequaloneDist_1=jequaloneDist_PT.(nm);
    end
    Params_ii=Params_solo{ii};
    Params_ii.r=PricePath.r(1);
    Params_ii.w=ParamPath.w(1);
    [~,Policy_ii]=ValueFnIter_Case1_FHorz(n_d_PT.(nm),n_a_PT.(nm),n_z_PT.(nm),N_j,d_grid_PT.(nm),a_grid_PT.(nm),z_grid_PT.(nm),pi_z_PT.(nm),ReturnFn_PT.(nm),Params_ii,DiscountFactorParamNames_GE,[],vfopt_solo{ii});
    AgentDist_initial_PT.(nm)=StationaryDist_FHorz_Case1(jequaloneDist_1,AgeWeightParamNames,Policy_ii,n_d_PT.(nm),n_a_PT.(nm),n_z_PT.(nm),N_j,pi_z_PT.(nm),Params_ii,simopt_solo{ii});
    Params_ii.r=PricePath.r(end);
    Params_ii.w=ParamPath.w(end);
    [V_final_PT.(nm),~]=ValueFnIter_Case1_FHorz(n_d_PT.(nm),n_a_PT.(nm),n_z_PT.(nm),N_j,d_grid_PT.(nm),a_grid_PT.(nm),z_grid_PT.(nm),pi_z_PT.(nm),ReturnFn_PT.(nm),Params_ii,DiscountFactorParamNames_GE,[],vfopt_solo{ii});
end
AgentDist_initial_PT.ptweights=[0.1;0.2;0.3;0.4];
PTypeDistParamNames={'ptweights_GEtest'};
Params_PT.ptweights_GEtest=[0.1;0.2;0.3;0.4];

%% One PType solve
transpathoptions.PricePathptype_vectoroutput=0;
[PricePath_GE,GEcondnPath_GE]=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist_PT, n_d_PT, n_a_PT, n_z_PT, N_j, Names_i, d_grid_PT,a_grid_PT,z_grid_PT, pi_z_PT, ReturnFn_PT, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames_GE, AgeWeightParamNames, PTypeDistParamNames, transpathoptions, simopt_PT, vfopt_PT);
transpathoptions.PricePathptype_vectoroutput=1;
PricePath_GEvec=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist_PT, n_d_PT, n_a_PT, n_z_PT, N_j, Names_i, d_grid_PT,a_grid_PT,z_grid_PT, pi_z_PT, ReturnFn_PT, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames_GE, AgeWeightParamNames, PTypeDistParamNames, transpathoptions, simopt_PT, vfopt_PT);
transpathoptions.PricePathptype_vectoroutput=0;

%% Four solo solves, without permanent types
transpathoptions_solo=transpathoptions;
transpathoptions_solo=rmfield(transpathoptions_solo,'GEptype');
transpathoptions_solo=rmfield(transpathoptions_solo,'PricePathptype_vectoroutput');
for ii=1:N_i
    nm=Names_i{ii};
    FnsToEvaluate_solo=struct();
    FnsToEvaluate_solo.assets=FnsToEvaluate.assets.(nm);
    FnsToEvaluate_solo.earnings=FnsToEvaluate.earnings.(nm);
    [PricePath_solo,GEcondnPath_solo]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final_PT.(nm), AgentDist_initial_PT.(nm), jequaloneDist_PT.(nm), n_d_PT.(nm), n_a_PT.(nm), n_z_PT.(nm), N_j, d_grid_PT.(nm),a_grid_PT.(nm),z_grid_PT.(nm), pi_z_PT.(nm), ReturnFn_PT.(nm), FnsToEvaluate_solo, GeneralEqmEqns, Params_solo{ii}, DiscountFactorParamNames_GE, AgeWeightParamNames, transpathoptions_solo, simopt_solo{ii}, vfopt_solo{ii});

    fprintf('GE decoupled types (fastOLG=%i), PType vs solo, PricePath r (type %s), this should be zero: %.3e \n',fastOLG,nm,max(abs(PricePath_GE.r.(nm)(:)-PricePath_solo.r(:))))
    fprintf('GE decoupled types (fastOLG=%i), PType vs solo, GEcondnPath (type %s), this should be zero: %.3e \n',fastOLG,nm,max(abs(GEcondnPath_GE.CapitalMarket.(nm)(:)-GEcondnPath_solo.CapitalMarket(:))))
    fprintf('GE decoupled types (fastOLG=%i), PricePathptype_vectoroutput 1 vs 0 (type %s), this should be zero: %.3e \n',fastOLG,nm,max(abs(reshape(PricePath_GEvec.r(ii,:),[],1)-PricePath_GE.r.(nm)(:))))
end

%% Anderson acceleration (GEnewprice=2)
% Anderson is not separable by ptype: each step solves one least-squares problem over the whole price
% path, every ptype's interest rate included, so it does not equal four solo Anderson solves. With
% anderson.memory=0 there is no history to mix and every step is a plain shooting step, so the PType
% Anderson solve must equal the PType shooting solve exactly. At maxiter=3 (iterations 2 and 3 would
% mix if there were a history, as warmup=0); and with the default memory, Anderson must then differ
% from shooting.
transpathoptionsA=rmfield(transpathoptions,{'GEnewprice3','PricePathptype_vectoroutput'});
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=transpathoptions.GEnewprice3.howtoupdate;
transpathoptionsA.maxiter=3;
transpathoptionsA.anderson.warmup=0;
transpathoptionsA.anderson.safeguard=0;
transpathoptionsA.anderson.memory=0;
PricePathA0=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist_PT, n_d_PT, n_a_PT, n_z_PT, N_j, Names_i, d_grid_PT,a_grid_PT,z_grid_PT, pi_z_PT, ReturnFn_PT, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames_GE, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsA, simopt_PT, vfopt_PT);
transpathoptionsA=rmfield(transpathoptionsA,'anderson');
transpathoptionsA.anderson.warmup=0;
transpathoptionsA.anderson.safeguard=0;
PricePathA5=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist_PT, n_d_PT, n_a_PT, n_z_PT, N_j, Names_i, d_grid_PT,a_grid_PT,z_grid_PT, pi_z_PT, ReturnFn_PT, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames_GE, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsA, simopt_PT, vfopt_PT);
transpathoptionsS3=rmfield(transpathoptions,'PricePathptype_vectoroutput');
transpathoptionsS3.maxiter=3;
PricePathS3=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist_PT, n_d_PT, n_a_PT, n_z_PT, N_j, Names_i, d_grid_PT,a_grid_PT,z_grid_PT, pi_z_PT, ReturnFn_PT, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames_GE, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsS3, simopt_PT, vfopt_PT);
maxdiff0=0; maxdiff5=0;
for ii=1:N_i
    nm=Names_i{ii};
    maxdiff0=max(maxdiff0,max(abs(PricePathA0.r.(nm)(:)-PricePathS3.r.(nm)(:))));
    maxdiff5=max(maxdiff5,max(abs(PricePathA5.r.(nm)(:)-PricePathS3.r.(nm)(:))));
end
fprintf('GE decoupled types (fastOLG=%i), Anderson with memory=0 vs shooting, both at maxiter=3, all four types, this should be zero: %.3e \n',fastOLG,maxdiff0)
fprintf('GE decoupled types (fastOLG=%i), Anderson with the default memory vs shooting, both at maxiter=3, all four types, this should NOT be zero: %.3e \n',fastOLG,maxdiff5)

output=struct();

end
