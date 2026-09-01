function output=EZRiskyAsset_CrossTests_semizasz_nod1_withA1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% CrossTest (semiz-as-z): RiskyAsset, nod1, WITH a1 -- EPSTEIN-ZIN version.
% a=[a1,a2]: a1=standard safe asset (a1prime chosen directly), a2=risky asset (set by aprimeFn).
% n_a=[n_a1,n_a2]; a_grid=[a1_grid;a2_grid].
% A semiz whose transitions do NOT depend on the decision (JustAMarkov) is just a markov.
% With searcheffortcost=0 (so dsemiz is inert) and uempbenefit=0, the semiz model reproduces
% the z-markov model exactly. V and Dist match bit-exact; Policy matches on [riskyshare,savings]
% (the inert dsemiz policy row is dropped, as its argmax is arbitrary).
% Mirrors CoreFHorzRiskyAsset_CrossTests_semizasz_nod1_withA1 mechanics, for the three EZ cases.
% TEST-FIRST: the semiz-side solves need EZ riskyasset semiz withA1 raws (do not exist yet) with
% n_z=0.
%
% Driver passes the SEMIZ n_d=[n_d2,n_d3,n_d4]=[riskyshare,savings,dsemiz]; d_grid=[d2;d3;d4].
% The z-side uses just the first two decisions.

% z-side decisions (drop the semiz decision dsemiz)
n_d_z=n_d(1:2);
d_grid_z=d_grid(1:(n_d(1)+n_d(2)));

% Make dsemiz inert and unemployment payoff zero, so semiz==z-markov
Params.searcheffortcost=0;
Params.uempbenefit=0;

% z-markov equal to the JustAMarkov semiz: grid [0;1], transitions from probfindjob/problosejob
z1=0; z2=1;
Params.z1=z1; Params.z2=z2;
z_grid=[z1;z2];
pi_z=[1-Params.probfindjob, Params.probfindjob; Params.problosejob, 1-Params.problosejob];
n_z=2;

SemiExoStateFn_JustAMarkov=@(n,nprime,dsemiz,probfindjob,problosejob,z1,z2) CoreFHorzRiskyAssetSetup_SemiExoStateFn_JustAMarkov(n,nprime,dsemiz,probfindjob,problosejob,z1,z2);

% ReturnFns: 4 shapes x 3 EZ cases
ReturnFn_z_noe_cons=@(savings,a1prime,a1,a2,z,r,w,kappa_j,r_a1,agej,Jr,pension) EZRiskyReturnFn_cons_nod1_z_noe_nosemiz_withA1(savings,a1prime,a1,a2,z,r,w,kappa_j,r_a1,agej,Jr,pension);
ReturnFn_s_noe_cons=@(savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_cons_nod1_noz_noe_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z_e_cons=@(savings,a1prime,a1,a2,z,e,r,w,kappa_j,r_a1,agej,Jr,pension) EZRiskyReturnFn_cons_nod1_z_e_nosemiz_withA1(savings,a1prime,a1,a2,z,e,r,w,kappa_j,r_a1,agej,Jr,pension);
ReturnFn_s_e_cons=@(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_cons_nod1_noz_e_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z_noe_posU=@(savings,a1prime,a1,a2,z,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension) EZRiskyReturnFn_positiveUtils_nod1_z_noe_nosemiz_withA1(savings,a1prime,a1,a2,z,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension);
ReturnFn_s_noe_posU=@(savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_positiveUtils_nod1_noz_noe_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z_e_posU=@(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension) EZRiskyReturnFn_positiveUtils_nod1_z_e_nosemiz_withA1(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension);
ReturnFn_s_e_posU=@(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_positiveUtils_nod1_noz_e_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z_noe_negU=@(savings,a1prime,a1,a2,z,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension) EZRiskyReturnFn_negativeUtils_nod1_z_noe_nosemiz_withA1(savings,a1prime,a1,a2,z,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension);
ReturnFn_s_noe_negU=@(savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_negativeUtils_nod1_noz_noe_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z_e_negU=@(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension) EZRiskyReturnFn_negativeUtils_nod1_z_e_nosemiz_withA1(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension);
ReturnFn_s_e_negU=@(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_negativeUtils_nod1_noz_e_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);

% Common riskyasset options (the EZ fields are added per case in the loop below)
vfoptions_common=struct();
vfoptions_common.riskyasset=1;
vfoptions_common.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions_common.n_u=vfoptionsbaseline.n_u; vfoptions_common.u_grid=vfoptionsbaseline.u_grid; vfoptions_common.pi_u=vfoptionsbaseline.pi_u;
simoptions=struct();
simoptions.riskyasset=1; simoptions.aprimeFn=vfoptions_common.aprimeFn;
simoptions.n_u=vfoptions_common.n_u; simoptions.u_grid=vfoptions_common.u_grid; simoptions.pi_u=vfoptions_common.pi_u;
simoptions.a_grid=a_grid;

% z-side options (2 decisions, refine_d=[0,1,1])
simoptions_z=simoptions; simoptions_z.refine_d=[0,1,1]; simoptions_z.d_grid=d_grid_z;
simoptions_z_e=simoptions_z; simoptions_z_e.n_e=simoptionsbaseline.n_e; simoptions_z_e.e_grid=simoptionsbaseline.e_grid; simoptions_z_e.pi_e=simoptionsbaseline.pi_e;

% semiz-side options (3 decisions, refine_d=[0,1,1,1], JustAMarkov semiz)
n_semiz=2;
simoptions_s=simoptions; simoptions_s.refine_d=[0,1,1,1]; simoptions_s.d_grid=d_grid;
simoptions_s.n_semiz=n_semiz; simoptions_s.semiz_grid=[z1;z2]; simoptions_s.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
simoptions_s_e=simoptions_s; simoptions_s_e.n_e=simoptionsbaseline.n_e; simoptions_s_e.e_grid=simoptionsbaseline.e_grid; simoptions_s_e.pi_e=simoptionsbaseline.pi_e;

%% The three EZ cases
for ezcase=1:3
    vfoptions=vfoptions_common;
    vfoptions.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptions.EZutils=0;
        vfoptions.EZriskaversion='ezgamma';
        vfoptions.EZeis='ezphi';
        ReturnFn_z_noe=ReturnFn_z_noe_cons; ReturnFn_s_noe=ReturnFn_s_noe_cons; ReturnFn_z_e=ReturnFn_z_e_cons; ReturnFn_s_e=ReturnFn_s_e_cons;
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptions.EZutils=1;
        vfoptions.EZpositiveutility=1;
        vfoptions.EZriskaversion='ezrisk';
        ReturnFn_z_noe=ReturnFn_z_noe_posU; ReturnFn_s_noe=ReturnFn_s_noe_posU; ReturnFn_z_e=ReturnFn_z_e_posU; ReturnFn_s_e=ReturnFn_s_e_posU;
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptions.EZutils=1;
        vfoptions.EZpositiveutility=0;
        vfoptions.EZriskaversion='ezrisk';
        ReturnFn_z_noe=ReturnFn_z_noe_negU; ReturnFn_s_noe=ReturnFn_s_noe_negU; ReturnFn_z_e=ReturnFn_z_e_negU; ReturnFn_s_e=ReturnFn_s_e_negU;
    end
    % z-side options (2 decisions, refine_d=[0,1,1])
    vfoptions_z=vfoptions; vfoptions_z.refine_d=[0,1,1];
    vfoptions_z_e=vfoptions_z; vfoptions_z_e.n_e=vfoptionsbaseline.n_e; vfoptions_z_e.e_grid=vfoptionsbaseline.e_grid; vfoptions_z_e.pi_e=vfoptionsbaseline.pi_e;
    % semiz-side options (3 decisions, refine_d=[0,1,1,1], JustAMarkov semiz)
    vfoptions_s=vfoptions; vfoptions_s.refine_d=[0,1,1,1];
    vfoptions_s.n_semiz=n_semiz; vfoptions_s.semiz_grid=[z1;z2]; vfoptions_s.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
    vfoptions_s_e=vfoptions_s; vfoptions_s_e.n_e=vfoptionsbaseline.n_e; vfoptions_s_e.e_grid=vfoptionsbaseline.e_grid; vfoptions_s_e.pi_e=vfoptionsbaseline.pi_e;

    %% (1) no e
    jequaloneDist_z=zeros([n_a,n_z],'gpuArray'); jequaloneDist_z(1,1,1)=1;
    jequaloneDist_s=zeros([n_a,n_semiz],'gpuArray'); jequaloneDist_s(1,1,1)=1;

    [V_z,Policy_z]=ValueFnIter_Case1_FHorz(n_d_z,n_a,n_z,N_j,d_grid_z,a_grid,z_grid,pi_z,ReturnFn_z_noe,Params,DiscountFactorParamNames,[],vfoptions_z);
    StationaryDist_z=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy_z,n_d_z,n_a,n_z,N_j,pi_z,Params,simoptions_z);

    [V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_s_noe,Params,DiscountFactorParamNames,[],vfoptions_s);
    StationaryDist_s=StationaryDist_FHorz_Case1(jequaloneDist_s,AgeWeightParamNames,Policy_s,n_d,n_a,0,N_j,[],Params,simoptions_s);

    Pd_z=Policy_z(1:2,:); Pd_s=Policy_s(1:2,:);
    fprintf('CrossTest semizasz (nod1 noe, withA1) [EZ %s]: should be zero: V %.3e, Policy(riskyshare,savings) %.3e, Dist %.3e \n',casestr,max(abs(V_z(:)-V_s(:))),max(abs(Pd_z(:)-Pd_s(:))),max(abs(StationaryDist_z(:)-StationaryDist_s(:))))

    %% (2) with e
    jequaloneDist_z=zeros([n_a,n_z,vfoptions_z_e.n_e],'gpuArray'); jequaloneDist_z(1,1,1,ceil(vfoptions_z_e.n_e/2))=1;
    jequaloneDist_s=zeros([n_a,n_semiz,vfoptions_s_e.n_e],'gpuArray'); jequaloneDist_s(1,1,1,ceil(vfoptions_s_e.n_e/2))=1;

    [V_z,Policy_z]=ValueFnIter_Case1_FHorz(n_d_z,n_a,n_z,N_j,d_grid_z,a_grid,z_grid,pi_z,ReturnFn_z_e,Params,DiscountFactorParamNames,[],vfoptions_z_e);
    StationaryDist_z=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy_z,n_d_z,n_a,n_z,N_j,pi_z,Params,simoptions_z_e);

    [V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_s_e,Params,DiscountFactorParamNames,[],vfoptions_s_e);
    StationaryDist_s=StationaryDist_FHorz_Case1(jequaloneDist_s,AgeWeightParamNames,Policy_s,n_d,n_a,0,N_j,[],Params,simoptions_s_e);

    Pd_z=Policy_z(1:2,:); Pd_s=Policy_s(1:2,:);
    fprintf('CrossTest semizasz (nod1 e, withA1) [EZ %s]: should be zero: V %.3e, Policy(riskyshare,savings) %.3e, Dist %.3e \n',casestr,max(abs(V_z(:)-V_s(:))),max(abs(Pd_z(:)-Pd_s(:))),max(abs(StationaryDist_z(:)-StationaryDist_s(:))))

end

%%
output=struct();

end
