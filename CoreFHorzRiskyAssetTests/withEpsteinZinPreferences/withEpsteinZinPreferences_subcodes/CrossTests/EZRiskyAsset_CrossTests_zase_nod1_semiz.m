function output=EZRiskyAsset_CrossTests_zase_nod1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% CrossTest (z-as-e) WITH a semiz background: RiskyAsset, nod1, single asset -- EPSTEIN-ZIN version.
% Same three z-as-e equivalences as EZRiskyAsset_CrossTests_zase_nod1, but every model also
% carries the (employment) semi-exogenous state. State-dimension ordering is (a, semiz, z, e).
% Driver passes the SEMIZ n_d=[riskyshare,savings,dsemiz]; refine_d=[0,1,1,1].
% Mirrors CoreFHorzRiskyAsset_CrossTests_zase_nod1_semiz mechanics, run for the three EZ cases.
% Under EZ the semiz transitions must also sit inside the ONE joint certainty-equivalent over
% (u,semizprime,zprime,eprime), so these tests are sharper than their vNM counterparts.
% TEST-FIRST: the legs that solve with n_z=0 ERROR until the EZ riskyasset semiz raws handle
% n_z=0 (and EZ riskyasset semiz beyond the single nod1 raw is itself test-first).

% z as a copy of e
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;

% ReturnFns: 4 shapes x 3 EZ cases
ReturnFn_none_cons=@(savings,dsemiz,a,semiz,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_cons_nod1_noz_noe_semiz(savings,dsemiz,a,semiz,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z_cons=@(savings,dsemiz,a,semiz,z,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_cons_nod1_z_noe_semiz(savings,dsemiz,a,semiz,z,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_e_cons=@(savings,dsemiz,a,semiz,e,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_cons_nod1_noz_e_semiz(savings,dsemiz,a,semiz,e,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_ze_cons=@(savings,dsemiz,a,semiz,z,e,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_cons_nod1_z_e_semiz(savings,dsemiz,a,semiz,z,e,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_none_posU=@(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_positiveUtils_nod1_noz_noe_semiz(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z_posU=@(savings,dsemiz,a,semiz,z,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_positiveUtils_nod1_z_noe_semiz(savings,dsemiz,a,semiz,z,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_e_posU=@(savings,dsemiz,a,semiz,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_positiveUtils_nod1_noz_e_semiz(savings,dsemiz,a,semiz,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_ze_posU=@(savings,dsemiz,a,semiz,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_positiveUtils_nod1_z_e_semiz(savings,dsemiz,a,semiz,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_none_negU=@(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_negativeUtils_nod1_noz_noe_semiz(savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z_negU=@(savings,dsemiz,a,semiz,z,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_negativeUtils_nod1_z_noe_semiz(savings,dsemiz,a,semiz,z,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_e_negU=@(savings,dsemiz,a,semiz,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_negativeUtils_nod1_noz_e_semiz(savings,dsemiz,a,semiz,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_ze_negU=@(savings,dsemiz,a,semiz,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_negativeUtils_nod1_z_e_semiz(savings,dsemiz,a,semiz,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);

% Riskyasset + semiz base options (the EZ fields are added per case in the loop below)
vfoptions_base=struct();
vfoptions_base.riskyasset=1;
vfoptions_base.refine_d=[0,1,1,1];
vfoptions_base.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions_base.n_u=vfoptionsbaseline.n_u; vfoptions_base.u_grid=vfoptionsbaseline.u_grid; vfoptions_base.pi_u=vfoptionsbaseline.pi_u;
vfoptions_base.n_semiz=vfoptionsbaseline.n_semiz; vfoptions_base.semiz_grid=vfoptionsbaseline.semiz_grid; vfoptions_base.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions=struct();
simoptions.riskyasset=1; simoptions.refine_d=vfoptions_base.refine_d; simoptions.aprimeFn=vfoptions_base.aprimeFn;
simoptions.n_u=vfoptions_base.n_u; simoptions.u_grid=vfoptions_base.u_grid; simoptions.pi_u=vfoptions_base.pi_u;
simoptions.n_semiz=vfoptions_base.n_semiz; simoptions.semiz_grid=vfoptions_base.semiz_grid; simoptions.SemiExoStateFn=vfoptions_base.SemiExoStateFn;
simoptions.d_grid=d_grid; simoptions.a_grid=a_grid;
n_semiz=vfoptions_base.n_semiz;
simoptions_withe=simoptions;
simoptions_withe.n_e=simoptionsbaseline.n_e; simoptions_withe.e_grid=simoptionsbaseline.e_grid; simoptions_withe.pi_e=simoptionsbaseline.pi_e;
simoptions_ze2=simoptions; simoptions_ze2.n_e=1; simoptions_ze2.e_grid=1; simoptions_ze2.pi_e=1;

%% The three EZ cases
for ezcase=1:3
    vfoptions=vfoptions_base;
    vfoptions.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptions.EZutils=0;
        vfoptions.EZriskaversion='ezgamma';
        vfoptions.EZeis='ezphi';
        ReturnFn_none=ReturnFn_none_cons; ReturnFn_z=ReturnFn_z_cons; ReturnFn_e=ReturnFn_e_cons; ReturnFn_ze=ReturnFn_ze_cons;
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptions.EZutils=1;
        vfoptions.EZpositiveutility=1;
        vfoptions.EZriskaversion='ezrisk';
        ReturnFn_none=ReturnFn_none_posU; ReturnFn_z=ReturnFn_z_posU; ReturnFn_e=ReturnFn_e_posU; ReturnFn_ze=ReturnFn_ze_posU;
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptions.EZutils=1;
        vfoptions.EZpositiveutility=0;
        vfoptions.EZriskaversion='ezrisk';
        ReturnFn_none=ReturnFn_none_negU; ReturnFn_z=ReturnFn_z_negU; ReturnFn_e=ReturnFn_e_negU; ReturnFn_ze=ReturnFn_ze_negU;
    end
    % Options carrying an e shock
    vfoptions_withe=vfoptions;
    vfoptions_withe.n_e=vfoptionsbaseline.n_e; vfoptions_withe.e_grid=vfoptionsbaseline.e_grid; vfoptions_withe.pi_e=vfoptionsbaseline.pi_e;
    vfoptions_ze2=vfoptions; vfoptions_ze2.n_e=1; vfoptions_ze2.e_grid=1; vfoptions_ze2.pi_e=1;

    %% (A) A single point z (value 1, prob 1) is the same as no z
    jequaloneDist_none=zeros([n_a,n_semiz],'gpuArray'); jequaloneDist_none(1,ceil(n_semiz/2))=1;

    [V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d,n_a,0,N_j,[],Params,simoptions);

    jequaloneDist_none_z=zeros([n_a,n_semiz,1],'gpuArray'); jequaloneDist_none_z(1,ceil(n_semiz/2),1)=1;
    [V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist0z=StationaryDist_FHorz_Case1(jequaloneDist_none_z,AgeWeightParamNames,Policy0z,n_d,n_a,1,N_j,1,Params,simoptions);

    fprintf('CrossTest zase+semiz (A) single-point z == no z [EZ %s], this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',casestr,max(abs(V0(:)-V0z(:))),max(abs(Policy0(:)-Policy0z(:))),max(abs(StationaryDist0(:)-StationaryDist0z(:))))

    %% (B) A markov-z that is really an iid == the same shock done as e
    jequaloneDist_z=zeros([n_a,n_semiz,n_z],'gpuArray'); jequaloneDist_z(1,ceil(n_semiz/2),ceil(n_z/2))=1;
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

    jequaloneDist_e=zeros([n_a,n_semiz,vfoptions_withe.n_e],'gpuArray'); jequaloneDist_e(1,ceil(n_semiz/2),ceil(vfoptions_withe.n_e/2))=1;
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_withe);
    StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist_e,AgeWeightParamNames,Policy2,n_d,n_a,0,N_j,[],Params,simoptions_withe);

    fprintf('CrossTest zase+semiz (B) iid-markov-z == e [EZ %s], this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',casestr,max(abs(V1(:)-V2(:))),max(abs(Policy1(:)-Policy2(:))),max(abs(StationaryDist1(:)-StationaryDist2(:))))

    %% (C) z&e code with the 'other' a single point reproduces the z-only model
    % z just a single point
    [V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_withe);
    jequaloneDist3=zeros([n_a,n_semiz,1,vfoptions_withe.n_e],'gpuArray'); jequaloneDist3(1,ceil(n_semiz/2),1,ceil(vfoptions_withe.n_e/2))=1;
    StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist3,AgeWeightParamNames,Policy3,n_d,n_a,1,N_j,1,Params,simoptions_withe);
    V3=squeeze(V3); Policy3=squeeze(Policy3); StationaryDist3=squeeze(StationaryDist3);

    fprintf('CrossTest zase+semiz (C) z&e with z=1 [EZ %s], this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',casestr,max(abs(V1(:)-V3(:))),max(abs(Policy1(:)-Policy3(:))),max(abs(StationaryDist1(:)-StationaryDist3(:))))

    % e just a single point
    [V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_ze2);
    jequaloneDist4=zeros([n_a,n_semiz,n_z,1],'gpuArray'); jequaloneDist4(1,ceil(n_semiz/2),ceil(n_z/2),1)=1;
    StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist4,AgeWeightParamNames,Policy4,n_d,n_a,n_z,N_j,pi_z,Params,simoptions_ze2);
    V4=squeeze(V4); Policy4=squeeze(Policy4); StationaryDist4=squeeze(StationaryDist4);

    fprintf('CrossTest zase+semiz (C) z&e with e=1 [EZ %s], this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',casestr,max(abs(V1(:)-V4(:))),max(abs(Policy1(:)-Policy4(:))),max(abs(StationaryDist1(:)-StationaryDist4(:))))

end

%%
output=struct();

end
