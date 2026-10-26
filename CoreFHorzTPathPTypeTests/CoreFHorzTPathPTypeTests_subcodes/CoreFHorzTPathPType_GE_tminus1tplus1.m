function output=CoreFHorzTPathPType_GE_tminus1tplus1(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,fastOLG)
% _tminus1 and _tplus1 on the general eqm transition path with permanent types. Two types (kappa_j
% differs, given as a structure), each with its own interest rate and a general eqm condition
% conditional on ptype (transpathoptions.GEptype), so the types do not interact and the PType solve
% must equal two solo solves without permanent types. The general eqm condition uses r_tplus1 and the
% aggregate assets_tminus1, and a FnsToEvaluate uses r_tminus1, r_tplus1 and w_tminus1 (so each ptype
% sees them on its own, inside its own path solve). The initial values depend on ptype, given as a
% structure (r) and as a vector of length N_i (assets); w is common. The interest rate by ptype is a
% structure with its fields in the opposite order to Names_i.
% Shooting for maxiter=2, with transpathoptions.fastOLG=fastOLG. Only z (no d, no e, no semiz).

Names_i={'low','high'};
N_i=2;

n_d=0;
d_grid=[];
vfoptions=struct();
simoptions=struct();

ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate.interest=@(aprime,a,z,r_tminus1,r_tplus1,w_tminus1) a*r_tminus1+aprime*r_tplus1+w_tminus1*z;

GeneralEqmEqns.CapitalMarket=@(r,r_tplus1,assets_tminus1,interest,earnings) r-0.02-0.005*(r_tplus1-r)-0.01*assets_tminus1/(1+earnings)+0.001*interest;

jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

Params_PT=Params;
Params_PT.kappa_j=struct();
Params_PT.kappa_j.low=Params.kappa_j;
Params_PT.kappa_j.high=1.4*Params.kappa_j;
PTypeDistParamNames={'ptweights_GEtest'};
Params_PT.ptweights_GEtest=[0.7;0.3];

initialvalues_r=[0.04,0.045];
initialvalues_assets=[0.5,0.6];

%% Per type: the initial agent distribution (stationary at the period-1 prices) and V_final (stationary at the final prices)
V_final_PT=struct(); Policy_final_PT=struct(); AgentDist_initial_PT=struct();
Params_solo=cell(N_i,1);
for ii=1:N_i
    Params_solo{ii}=Params;
    Params_solo{ii}.kappa_j=Params_PT.kappa_j.(Names_i{ii});
    Params_ii=Params_solo{ii};
    Params_ii.r=PricePath.r(1);
    Params_ii.w=ParamPath.w(1);
    [~,Policy_ii]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initial_PT.(Names_i{ii})=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_ii,n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions);
    Params_ii.r=PricePath.r(end);
    Params_ii.w=ParamPath.w(end);
    [V_final_PT.(Names_i{ii}),Policy_final_PT.(Names_i{ii})]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
end
AgentDist_initial_PT.ptweights=[0.7;0.3];

%% One PType solve
% A price that depends on ptype, as a structure, deliberately with its fields in the opposite order to
% Names_i: the solver must read them by name
PricePath_PT=struct();
PricePath_PT.r=struct('high',PricePath.r,'low',PricePath.r);

transpathoptions=struct();
transpathoptions.fastOLG=fastOLG;
transpathoptions.GEnewprice=3; % the shooting algorithm
transpathoptions.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5};
transpathoptions.GEptype={'CapitalMarket'};
transpathoptions.maxiter=2; % returns the price path after one update, and the general eqm conditions evaluated on it
transpathoptions.verbose=0;
transpathoptions.initialvalues.r.low=initialvalues_r(1);
transpathoptions.initialvalues.r.high=initialvalues_r(2);
transpathoptions.initialvalues.w=1.2;
transpathoptions.initialvalues.assets=initialvalues_assets'; % one per ptype

[PricePath_GE,GEcondnPath_GE]=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);

%% Two solo solves, without permanent types
for ii=1:N_i
    nm=Names_i{ii};
    transpathoptions_solo=struct();
    transpathoptions_solo.fastOLG=fastOLG;
    transpathoptions_solo.GEnewprice=3;
    transpathoptions_solo.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5};
    transpathoptions_solo.maxiter=2;
    transpathoptions_solo.verbose=0;
    transpathoptions_solo.initialvalues.r=initialvalues_r(ii);
    transpathoptions_solo.initialvalues.w=1.2;
    transpathoptions_solo.initialvalues.assets=initialvalues_assets(ii);
    [PricePath_solo,GEcondnPath_solo]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final_PT.(nm), AgentDist_initial_PT.(nm), jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_solo{ii}, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions_solo, simoptions, vfoptions);

    fprintf('GE _tminus1/_tplus1 by ptype (fastOLG=%i), PType vs solo, PricePath r (type %s), this should be zero: %.3e \n',fastOLG,nm,max(abs(PricePath_GE.r.(nm)(:)-PricePath_solo.r(:))))
    fprintf('GE _tminus1/_tplus1 by ptype (fastOLG=%i), PType vs solo, GEcondnPath (type %s), this should be zero: %.3e \n',fastOLG,nm,max(abs(GEcondnPath_GE.CapitalMarket.(nm)(:)-GEcondnPath_solo.CapitalMarket(:))))
end

%% Anderson acceleration (GEnewprice=2)
% Anderson is not separable by ptype (each step is one least-squares problem over the whole price path),
% so it is not compared with solo solves. With anderson.memory=0 every step is a plain shooting step, so
% it must equal shooting exactly, at maxiter=3; with the default memory (and warmup=0) it must differ.
transpathoptionsA=rmfield(transpathoptions,'GEnewprice3');
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=transpathoptions.GEnewprice3.howtoupdate;
transpathoptionsA.maxiter=3;
transpathoptionsA.anderson.warmup=0;
transpathoptionsA.anderson.safeguard=0;
transpathoptionsA.anderson.memory=0;
PricePathA0=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsA, simoptions, vfoptions);
transpathoptionsA=rmfield(transpathoptionsA,'anderson');
transpathoptionsA.anderson.warmup=0;
transpathoptionsA.anderson.safeguard=0;
PricePathA5=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsA, simoptions, vfoptions);
transpathoptionsS3=transpathoptions;
transpathoptionsS3.maxiter=3;
PricePathS3=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsS3, simoptions, vfoptions);
fprintf('GE _tminus1/_tplus1 by ptype (fastOLG=%i), Anderson with memory=0 vs shooting, both at maxiter=3, both types, this should be zero: %.3e \n',fastOLG,max([max(abs(PricePathA0.r.low(:)-PricePathS3.r.low(:))),max(abs(PricePathA0.r.high(:)-PricePathS3.r.high(:)))]))
fprintf('GE _tminus1/_tplus1 by ptype (fastOLG=%i), Anderson with the default memory vs shooting, both at maxiter=3, both types, this should NOT be zero: %.3e \n',fastOLG,max([max(abs(PricePathA5.r.low(:)-PricePathS3.r.low(:))),max(abs(PricePathA5.r.high(:)-PricePathS3.r.high(:)))]))

%% (PA6) Statistics along the path with _tminus1 and _tplus1 by ptype
% Only for fastOLG=1 (the stats commands take the paths, laid out the same whatever the fastOLG).
% The standalone PType paths at the input price path, and AllStats and the life-cycle profiles with
% FnsToEvaluate that use r_tminus1, r_tplus1, w_tminus1 (prices and a parameter) and assets_tminus1 (an
% aggregate, which is the mean conditional on the ptype), with the initial values by ptype as above.
% Each ptype must equal the commands without permanent types on its own paths, with its own initial values.
if fastOLG==1
    OutputTypes={'Mean','Median','Variance','LorenzCurve','Gini','QuantileCutoffs','QuantileMeans'};
    FnsToEvaluateL.assets=FnsToEvaluate.assets;
    FnsToEvaluateL.interest=FnsToEvaluate.interest;
    FnsToEvaluateL.y=@(aprime,a,z,assets_tminus1) a*assets_tminus1;
    transpathoptions_sa=struct();
    transpathoptions_sa.fastOLG=1;
    transpathoptions_sa.initialvalues=transpathoptions.initialvalues;
    [~,PolicyPath_PT]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, Policy_final_PT, Params_PT, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
    AgentDistPath_PT=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist, PricePath_PT, ParamPath, PolicyPath_PT, AgeWeightParamNames, PTypeDistParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z, T,Params_PT, transpathoptions_sa, simoptions);
    AllStats_PT=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsToEvaluateL, AgentDistPath_PT, PolicyPath_PT, PricePath_PT, ParamPath, Params_PT, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);
    LCP_PT=LifeCycleProfiles_TransPath_FHorz_Case1_PType(FnsToEvaluateL, AgentDistPath_PT, PolicyPath_PT, PricePath_PT, ParamPath, Params_PT, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);
    fnnames=fieldnames(FnsToEvaluateL);
    for ii=1:N_i
        nm=Names_i{ii};
        transpathoptions_solo=struct();
        transpathoptions_solo.fastOLG=1;
        transpathoptions_solo.initialvalues.r=initialvalues_r(ii);
        transpathoptions_solo.initialvalues.w=1.2;
        transpathoptions_solo.initialvalues.assets=initialvalues_assets(ii);
        REF=EvalFnOnTransPath_AllStats_Case1_FHorz(FnsToEvaluateL, AgentDistPath_PT.(nm), PolicyPath_PT.(nm), PricePath, ParamPath, Params_solo{ii}, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptions_solo, simoptions);
        REFL=LifeCycleProfiles_TransPath_FHorz_Case1(FnsToEvaluateL, AgentDistPath_PT.(nm), PolicyPath_PT.(nm), PricePath, ParamPath, Params_solo{ii}, T, n_d, n_a, n_z, N_j, d_grid, a_grid, z_grid, transpathoptions_solo, simoptions);
        % every field of the reference, to any depth, found in the PType output with the ptype name after the FnsToEvaluate name
        maxdiff=0; nanmismatch=0;
        queue={};
        REFnames=fieldnames(REF);
        for kk=1:length(REFnames)
            queue{end+1}=REFnames(kk);
        end
        while ~isempty(queue)
            fieldpath=queue{1};
            queue(1)=[];
            O=getfield(REF,fieldpath{:});
            if isstruct(O)
                Onames=fieldnames(O);
                for kk=1:length(Onames)
                    queue{end+1}=[fieldpath,Onames(kk)];
                end
            else
                Ppath=[fieldpath(1),{nm},fieldpath(2:end)];
                P=AllStats_PT;
                found=1;
                for kk=1:length(Ppath)
                    if isstruct(P) && isfield(P,Ppath{kk})
                        P=P.(Ppath{kk});
                    else
                        found=0;
                        break
                    end
                end
                if found==0
                    nanmismatch=nanmismatch+1;
                elseif iscell(O) % text (the comments): must be the same text
                    nanmismatch=nanmismatch+(~isequal(P,O));
                else
                    maxdiff=max([maxdiff; abs(P(:)-O(:))],[],'omitnan');
                    nanmismatch=nanmismatch+sum(isnan(P(:))~=isnan(O(:)));
                end
            end
        end
        fprintf('PType stats, _tminus1/_tplus1 by ptype: (PA6) %s, EvalFnOnTransPath_AllStats_PType vs without permanent types with its own initial values, every stat, this should be zero: %.3e \n',nm,maxdiff)
        fprintf('PType stats, _tminus1/_tplus1 by ptype: (PA6) %s, EvalFnOnTransPath_AllStats_PType vs without permanent types with its own initial values, positions NaN (or text) in one but not the other, or missing, this should be zero: %i \n',nm,nanmismatch)
        maxdiff=0; nanmismatch=0;
        for ff=1:length(fnnames)
            for oo=1:length(OutputTypes)
                PI=gather(LCP_PT.(fnnames{ff}).(nm).initialages.(OutputTypes{oo}));
                PB=gather(LCP_PT.(fnnames{ff}).(nm).bornduringtranstion.(OutputTypes{oo}));
                OI=gather(REFL.(fnnames{ff}).initialages.(OutputTypes{oo}));
                OB=gather(REFL.(fnnames{ff}).bornduringtranstion.(OutputTypes{oo}));
                maxdiff=max([maxdiff; abs(PI(:)-OI(:)); abs(PB(:)-OB(:))],[],'omitnan');
                nanmismatch=nanmismatch+sum(isnan(PI(:))~=isnan(OI(:)))+sum(isnan(PB(:))~=isnan(OB(:)));
            end
        end
        fprintf('PType stats, _tminus1/_tplus1 by ptype: (PA6) %s, LifeCycleProfiles_TransPath_PType vs without permanent types with its own initial values, all seven outputs, this should be zero: %.3e \n',nm,maxdiff)
        fprintf('PType stats, _tminus1/_tplus1 by ptype: (PA6) %s, LifeCycleProfiles_TransPath_PType vs without permanent types with its own initial values, entries NaN in one but not the other, this should be zero: %i \n',nm,nanmismatch)
    end
    % (the initial values differ by ptype, so a command that gave every ptype the same one would fail the above)
end

output=struct();

end
