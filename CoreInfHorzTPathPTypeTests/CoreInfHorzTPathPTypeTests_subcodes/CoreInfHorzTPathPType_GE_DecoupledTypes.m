function output=CoreInfHorzTPathPType_GE_DecoupledTypes(T,PricePath,ParamPath,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline)
% The general eqm transition path with permanent types, for four types that differ in nearly everything
% a type can differ in, and whose capital market condition is conditional on ptype
% (transpathoptions.GEptype), with an interest rate per ptype. The types then do not interact, so the one
% PType solve must give exactly what four separate solves without permanent types give.
%   none: no z, no e; its own n_a and a_grid
%   z:    z; divide-and-conquer
%   e:    e; a decision variable d
%   ze:   z and e, with a z_grid different from type z (same n_z, other values)
% There is also a general eqm condition that is common to all ptypes, TaxBalance, in a price tau that no
% ptype uses, listed AFTER the one by ptype. Every solo solve has it too. With the by-ptype condition
% taking N_i columns of the internal GEcondnPath and the common one a single column, a GEcondnPath that
% is read back by position rather than by name puts the wrong numbers under the wrong names.
% The interest rate by ptype is a structure with its fields in the opposite order to Names_i.
% Each solve runs the shooting algorithm for maxiter=2 (so returns the price path after one update, and
% the general eqm conditions on it). Checked: PricePath and GEcondnPath of each type against its solo
% solve, the common condition and price, PricePathptype_vectoroutput=1 against =0, and Anderson with
% memory=0 against shooting (Anderson is not separable by ptype, so not against the solo solves).

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
ReturnFn_PT=struct();
vfopt_PT=struct(); simopt_PT=struct();
vfopt_PT.n_e=struct(); vfopt_PT.e_grid=struct(); vfopt_PT.pi_e=struct();
simopt_PT.n_e=struct(); simopt_PT.e_grid=struct(); simopt_PT.pi_e=struct();
vfopt_PT.divideandconquer=struct();

% type 1: none
nm='none';
n_d_PT.(nm)=0; d_grid_PT.(nm)=[];
n_a_PT.(nm)=n_a_none; a_grid_PT.(nm)=a_grid_none;
n_z_PT.(nm)=0; z_grid_PT.(nm)=[]; pi_z_PT.(nm)=[];
ReturnFn_PT.(nm)=@(aprime,a,r,w,sigma) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,sigma);
FnsToEvaluate.assets.(nm)=@(aprime,a) a;
FnsToEvaluate.earnings.(nm)=@(aprime,a,w) w;

% type 2: z, divide-and-conquer
nm='z';
n_d_PT.(nm)=0; d_grid_PT.(nm)=[];
n_a_PT.(nm)=n_a; a_grid_PT.(nm)=a_grid;
n_z_PT.(nm)=n_z; z_grid_PT.(nm)=z_grid; pi_z_PT.(nm)=pi_z;
ReturnFn_PT.(nm)=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
FnsToEvaluate.assets.(nm)=@(aprime,a,z) a;
FnsToEvaluate.earnings.(nm)=@(aprime,a,z,w) w*z;
vfopt_PT.divideandconquer.(nm)=1;

% type 3: e, with d
nm='e';
n_d_PT.(nm)=n_d; d_grid_PT.(nm)=d_grid;
n_a_PT.(nm)=n_a; a_grid_PT.(nm)=a_grid;
n_z_PT.(nm)=0; z_grid_PT.(nm)=[]; pi_z_PT.(nm)=[];
ReturnFn_PT.(nm)=@(d,aprime,a,e,r,w,sigma,eta,varphi) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,sigma,eta,varphi);
FnsToEvaluate.assets.(nm)=@(d,aprime,a,e) a;
FnsToEvaluate.earnings.(nm)=@(d,aprime,a,e,w) w*d*e;
vfopt_PT.n_e.(nm)=n_e; vfopt_PT.e_grid.(nm)=e_grid; vfopt_PT.pi_e.(nm)=pi_e;
simopt_PT.n_e.(nm)=n_e; simopt_PT.e_grid.(nm)=e_grid; simopt_PT.pi_e.(nm)=pi_e;

% type 4: z and e, z_grid differs from type z
nm='ze';
n_d_PT.(nm)=0; d_grid_PT.(nm)=[];
n_a_PT.(nm)=n_a; a_grid_PT.(nm)=a_grid;
n_z_PT.(nm)=n_z; z_grid_PT.(nm)=z_grid.^1.5; pi_z_PT.(nm)=pi_z;
ReturnFn_PT.(nm)=@(aprime,a,z,e,r,w,sigma) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,sigma);
FnsToEvaluate.assets.(nm)=@(aprime,a,z,e) a;
FnsToEvaluate.earnings.(nm)=@(aprime,a,z,e,w) w*z*e;
vfopt_PT.n_e.(nm)=n_e; vfopt_PT.e_grid.(nm)=e_grid; vfopt_PT.pi_e.(nm)=pi_e;
simopt_PT.n_e.(nm)=n_e; simopt_PT.e_grid.(nm)=e_grid; simopt_PT.pi_e.(nm)=pi_e;

%% General eqm: one condition per ptype, each with its own interest rate, then one common condition
GeneralEqmEqns.CapitalMarket=@(r,assets,earnings) r-0.02-0.01*assets/(1+earnings);
GeneralEqmEqns.TaxBalance=@(tau) tau-0.1;
Params.tau=0.12;
% The interest rate depends on ptype: a structure, with its fields in the opposite order to Names_i
PricePath_PT=struct();
PricePath_PT.r=struct('ze',PricePath.r,'e',PricePath.r,'z',PricePath.r,'none',PricePath.r);
PricePath_PT.tau=0.12*ones(1,T);
PricePath_solo=struct();
PricePath_solo.r=PricePath.r;
PricePath_solo.tau=PricePath_PT.tau;

transpathoptions=struct();
transpathoptions.GEnewprice=3; % the shooting algorithm
transpathoptions.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5; 'TaxBalance','tau',0,0.5};
transpathoptions.GEptype={'CapitalMarket'};
transpathoptions.maxiter=2; % returns the price path after one update, and the general eqm conditions evaluated on it
transpathoptions.verbose=0;

%% Each type: the initial agent distribution (stationary at the period-1 prices) and V_final (stationary at the final prices)
V_final_PT=struct(); AgentDist_initial_PT=struct();
vfopt_solo=cell(N_i,1); simopt_solo=cell(N_i,1);
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
    Params_ii=Params;
    Params_ii.r=PricePath.r(1);
    Params_ii.w=ParamPath.w(1);
    [~,Policy_ii]=ValueFnIter_InfHorz(n_d_PT.(nm),n_a_PT.(nm),n_z_PT.(nm),d_grid_PT.(nm),a_grid_PT.(nm),z_grid_PT.(nm),pi_z_PT.(nm),ReturnFn_PT.(nm),Params_ii,DiscountFactorParamNames,[],vfopt_solo{ii});
    AgentDist_initial_PT.(nm)=StationaryDist_InfHorz(Policy_ii,n_d_PT.(nm),n_a_PT.(nm),n_z_PT.(nm),pi_z_PT.(nm),simopt_solo{ii},Params_ii,[]);
    Params_ii.r=PricePath.r(end);
    Params_ii.w=ParamPath.w(end);
    [V_final_PT.(nm),~]=ValueFnIter_InfHorz(n_d_PT.(nm),n_a_PT.(nm),n_z_PT.(nm),d_grid_PT.(nm),a_grid_PT.(nm),z_grid_PT.(nm),pi_z_PT.(nm),ReturnFn_PT.(nm),Params_ii,DiscountFactorParamNames,[],vfopt_solo{ii});
end
AgentDist_initial_PT.ptweights=[0.1;0.2;0.3;0.4];
PTypeDistParamNames={'ptweights_GEtest'};
Params_PT=Params;
Params_PT.ptweights_GEtest=[0.1;0.2;0.3;0.4];

%% One PType solve
transpathoptions.PricePathptype_vectoroutput=0;
[PricePath_GE,GEcondnPath_GE]=TransitionPath_InfHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, n_d_PT,n_a_PT,n_z_PT, Names_i, d_grid_PT,a_grid_PT,z_grid_PT, pi_z_PT, ReturnFn_PT, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions, simopt_PT, vfopt_PT);
transpathoptions.PricePathptype_vectoroutput=1;
PricePath_GEvec=TransitionPath_InfHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, n_d_PT,n_a_PT,n_z_PT, Names_i, d_grid_PT,a_grid_PT,z_grid_PT, pi_z_PT, ReturnFn_PT, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions, simopt_PT, vfopt_PT);
transpathoptions.PricePathptype_vectoroutput=0;

%% Four solo solves, without permanent types
transpathoptions_solo=rmfield(transpathoptions,{'GEptype','PricePathptype_vectoroutput'});
for ii=1:N_i
    nm=Names_i{ii};
    FnsToEvaluate_solo=struct();
    FnsToEvaluate_solo.assets=FnsToEvaluate.assets.(nm);
    FnsToEvaluate_solo.earnings=FnsToEvaluate.earnings.(nm);
    [PricePath_s,GEcondnPath_s]=TransitionPath_InfHorz(PricePath_solo, ParamPath, T, V_final_PT.(nm), AgentDist_initial_PT.(nm), n_d_PT.(nm),n_a_PT.(nm),n_z_PT.(nm), d_grid_PT.(nm),a_grid_PT.(nm),z_grid_PT.(nm), pi_z_PT.(nm), ReturnFn_PT.(nm), FnsToEvaluate_solo, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions_solo, simopt_solo{ii}, vfopt_solo{ii}, []);
    GEs=gather(GEcondnPath_s.CapitalMarket(:));
    GEp=gather(GEcondnPath_GE.CapitalMarket.(nm)(:));
    fprintf('GE decoupled types, PType vs solo, PricePath r (type %s), this should be zero: %.3e \n',nm,max(abs(gather(PricePath_GE.r.(nm)(:))-gather(PricePath_s.r(:)))))
    fprintf('GE decoupled types, PType vs solo, GEcondnPath CapitalMarket (type %s), this should be zero: %.3e \n',nm,max(abs(GEp(1:T-1)-GEs(1:T-1))))
    fprintf('GE decoupled types, PricePathptype_vectoroutput 1 vs 0 (type %s), this should be zero: %.3e \n',nm,max(abs(gather(reshape(PricePath_GEvec.r(ii,:),[],1))-gather(PricePath_GE.r.(nm)(:)))))
    if ii==1
        GEts=gather(GEcondnPath_s.TaxBalance(:));
        GEtp=gather(GEcondnPath_GE.TaxBalance(:));
        fprintf('GE decoupled types, PType vs solo, common condition, GEcondnPath TaxBalance, this should be zero: %.3e \n',max(abs(GEtp(1:T-1)-GEts(1:T-1))))
        fprintf('GE decoupled types, PType vs solo, common price, PricePath tau, this should be zero: %.3e \n',max(abs(gather(PricePath_GE.tau(:))-gather(PricePath_s.tau(:)))))
    end
end
fprintf('GE decoupled types, for scale the update moved r of type ze by %.3e (must not be zero, or the PricePath checks above are empty) \n',max(abs(gather(PricePath_GE.r.ze(:))-PricePath.r(:))))

%% Anderson acceleration (GEnewprice=2)
% Anderson is not separable by ptype: each step solves one least-squares problem over the whole price
% path, every ptype's interest rate included, so it does not equal four solo Anderson solves. With
% anderson.memory=0 there is no history to mix and every step is a plain shooting step, so the PType
% Anderson solve must equal the PType shooting solve exactly, at maxiter=3 (iterations 2 and 3 would
% mix if there were a history, as warmup=0). With the default memory it must then differ from shooting.
transpathoptionsA=rmfield(transpathoptions,{'GEnewprice3','PricePathptype_vectoroutput'});
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=transpathoptions.GEnewprice3.howtoupdate;
transpathoptionsA.maxiter=3;
transpathoptionsA.anderson.warmup=0;
transpathoptionsA.anderson.safeguard=0;
transpathoptionsA.anderson.memory=0;
PricePathA0=TransitionPath_InfHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, n_d_PT,n_a_PT,n_z_PT, Names_i, d_grid_PT,a_grid_PT,z_grid_PT, pi_z_PT, ReturnFn_PT, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsA, simopt_PT, vfopt_PT);
transpathoptionsA=rmfield(transpathoptionsA,'anderson');
transpathoptionsA.anderson.warmup=0;
transpathoptionsA.anderson.safeguard=0;
PricePathA5=TransitionPath_InfHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, n_d_PT,n_a_PT,n_z_PT, Names_i, d_grid_PT,a_grid_PT,z_grid_PT, pi_z_PT, ReturnFn_PT, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsA, simopt_PT, vfopt_PT);
transpathoptionsS3=rmfield(transpathoptions,'PricePathptype_vectoroutput');
transpathoptionsS3.maxiter=3;
PricePathS3=TransitionPath_InfHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, n_d_PT,n_a_PT,n_z_PT, Names_i, d_grid_PT,a_grid_PT,z_grid_PT, pi_z_PT, ReturnFn_PT, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsS3, simopt_PT, vfopt_PT);
maxdiff0=max(abs(gather(PricePathA0.tau(:))-gather(PricePathS3.tau(:))));
maxdiff5=max(abs(gather(PricePathA5.tau(:))-gather(PricePathS3.tau(:))));
for ii=1:N_i
    nm=Names_i{ii};
    maxdiff0=max(maxdiff0,max(abs(gather(PricePathA0.r.(nm)(:))-gather(PricePathS3.r.(nm)(:)))));
    maxdiff5=max(maxdiff5,max(abs(gather(PricePathA5.r.(nm)(:))-gather(PricePathS3.r.(nm)(:)))));
end
fprintf('GE decoupled types, Anderson with memory=0 vs shooting, both at maxiter=3, all four types and tau, this should be zero: %.3e \n',maxdiff0)
fprintf('GE decoupled types, Anderson with the default memory vs shooting, both at maxiter=3, all four types and tau, this should NOT be zero: %.3e \n',maxdiff5)

output=struct();

end
