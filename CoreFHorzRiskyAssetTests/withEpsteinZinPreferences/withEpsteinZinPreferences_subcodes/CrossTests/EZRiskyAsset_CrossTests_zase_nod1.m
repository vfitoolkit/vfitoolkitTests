function output=EZRiskyAsset_CrossTests_zase_nod1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% CrossTest (z-as-e): RiskyAsset, nod1, nosemiz, single asset -- EPSTEIN-ZIN version.
% Mirrors CoreFHorzRiskyAsset_CrossTests_zase_nod1 mechanics, run for the three EZ cases
% (cons-units / positive utils / negative utils), following the main EZ bank convention.
% Three internal equivalences (all riskyasset, all should be machine-precision zero):
%   (A) a single-point z (value 1, prob 1) gives the same as no z at all
%   (B) a markov-z that is really an iid gives the same as the same shock done as e (iid)
%   (C) solving the z&e code with the 'other' shock a single point reproduces the z-only model
% Under EZ these are NOT redundant re-checks of vNM machinery: with the u shock present they are
% the riskyasset analogue of EZFHorz_CrossTests3 (nested-vs-joint certainty-equivalents). The
% same shock must enter the SAME single joint CE over (u,zprime,eprime) whether declared as z or
% as e; nested CEs == joint CE only at the gamma=1/phi knife-edge, and the defaults (ezgamma=3,
% ezphi=0.5, ezrisk=3) are away from it.
% TEST-FIRST: the legs that solve with n_z=0 (block (A) no-z side, block (B) e-side, block (C)
% z=1 side runs n_z=1 so is fine) ERROR until the EZ riskyasset raws handle n_z=0.

% For crosstests, set z to be a copy of e
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;
% NOTE: z & e appear in the same place in the ReturnFn

% ReturnFns: 4 shapes x 3 EZ cases
ReturnFn_none_cons=@(savings,a,r,w,kappa_j,agej,Jr,pension) EZRiskyReturnFn_cons_nod1_noz_noe_nosemiz(savings,a,r,w,kappa_j,agej,Jr,pension);
ReturnFn_z_cons=@(savings,a,z,r,w,kappa_j,agej,Jr,pension) EZRiskyReturnFn_cons_nod1_z_noe_nosemiz(savings,a,z,r,w,kappa_j,agej,Jr,pension);
ReturnFn_e_cons=@(savings,a,e,r,w,kappa_j,agej,Jr,pension) EZRiskyReturnFn_cons_nod1_noz_e_nosemiz(savings,a,e,r,w,kappa_j,agej,Jr,pension);
ReturnFn_ze_cons=@(savings,a,z,e,r,w,kappa_j,agej,Jr,pension) EZRiskyReturnFn_cons_nod1_z_e_nosemiz(savings,a,z,e,r,w,kappa_j,agej,Jr,pension);
ReturnFn_none_posU=@(savings,a,r,w,kappa_j,ezsigma,agej,Jr,pension) EZRiskyReturnFn_positiveUtils_nod1_noz_noe_nosemiz(savings,a,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_z_posU=@(savings,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension) EZRiskyReturnFn_positiveUtils_nod1_z_noe_nosemiz(savings,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_e_posU=@(savings,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension) EZRiskyReturnFn_positiveUtils_nod1_noz_e_nosemiz(savings,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_ze_posU=@(savings,a,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension) EZRiskyReturnFn_positiveUtils_nod1_z_e_nosemiz(savings,a,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_none_negU=@(savings,a,r,w,kappa_j,ezsigma,agej,Jr,pension) EZRiskyReturnFn_negativeUtils_nod1_noz_noe_nosemiz(savings,a,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_z_negU=@(savings,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension) EZRiskyReturnFn_negativeUtils_nod1_z_noe_nosemiz(savings,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_e_negU=@(savings,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension) EZRiskyReturnFn_negativeUtils_nod1_noz_e_nosemiz(savings,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_ze_negU=@(savings,a,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension) EZRiskyReturnFn_negativeUtils_nod1_z_e_nosemiz(savings,a,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension);

% Riskyasset base options (the EZ fields are added per case in the loop below)
vfoptions_base=struct();
vfoptions_base.riskyasset=1;
vfoptions_base.refine_d=[0,1,1];
vfoptions_base.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions_base.n_u=vfoptionsbaseline.n_u; vfoptions_base.u_grid=vfoptionsbaseline.u_grid; vfoptions_base.pi_u=vfoptionsbaseline.pi_u;
simoptions=struct();
simoptions.riskyasset=1; simoptions.refine_d=vfoptions_base.refine_d; simoptions.aprimeFn=vfoptions_base.aprimeFn;
simoptions.n_u=vfoptions_base.n_u; simoptions.u_grid=vfoptions_base.u_grid; simoptions.pi_u=vfoptions_base.pi_u;
simoptions.d_grid=d_grid; simoptions.a_grid=a_grid;
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

    %% (A) A single point z (value 1, prob 1) is the same as no shocks
    jequaloneDist_none=zeros([n_a,1],'gpuArray'); jequaloneDist_none(1,1,1)=1;

    [V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d,n_a,0,N_j,[],Params,simoptions);

    [V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist0z=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0z,n_d,n_a,1,N_j,1,Params,simoptions);

    fprintf('CrossTest zase (A) single-point z == no z [EZ %s], this should be zero: V %.3e, Policy %.3e, Dist %.3e \n',casestr,max(abs(V0(:)-V0z(:))),max(abs(Policy0(:)-Policy0z(:))),max(abs(StationaryDist0(:)-StationaryDist0z(:))))

    %% (B) A markov-z that is really an iid == the same shock done as e
    jequaloneDist_z=zeros([n_a,n_z],'gpuArray'); jequaloneDist_z(1,ceil(n_z/2))=1;

    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_withe);
    StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy2,n_d,n_a,0,N_j,[],Params,simoptions_withe);

    fprintf('CrossTest zase (B) iid-markov-z == e [EZ %s], this should be zero: V %.3e, Policy %.3e, Dist %.3e \n',casestr,max(abs(V1(:)-V2(:))),max(abs(Policy1(:)-Policy2(:))),max(abs(StationaryDist1(:)-StationaryDist2(:))))

    %% (C) z&e code with the 'other' a single point reproduces the z-only model
    % First, make z just a single point (so only e is active, but e is a copy of z)
    [V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_withe);
    jequaloneDist3=zeros([n_a,1,vfoptions_withe.n_e],'gpuArray'); jequaloneDist3(1,1,ceil(vfoptions_withe.n_e/2))=1;
    StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist3,AgeWeightParamNames,Policy3,n_d,n_a,1,N_j,1,Params,simoptions_withe);
    V3=squeeze(V3); Policy3=squeeze(Policy3); StationaryDist3=squeeze(StationaryDist3);

    fprintf('CrossTest zase (C) z&e with z=1 [EZ %s], this should be zero: V %.3e, Policy %.3e, Dist %.3e \n',casestr,max(abs(V1(:)-V3(:))),max(abs(Policy1(:)-Policy3(:))),max(abs(StationaryDist1(:)-StationaryDist3(:))))

    % Second, make e just a single point (so only z is active)
    [V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_ze2);
    jequaloneDist4=zeros([n_a,n_z,1],'gpuArray'); jequaloneDist4(1,ceil(n_z/2),1)=1;
    StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist4,AgeWeightParamNames,Policy4,n_d,n_a,n_z,N_j,pi_z,Params,simoptions_ze2);
    V4=squeeze(V4); Policy4=squeeze(Policy4); StationaryDist4=squeeze(StationaryDist4);

    fprintf('CrossTest zase (C) z&e with e=1 [EZ %s], this should be zero: V %.3e, Policy %.3e, Dist %.3e \n',casestr,max(abs(V1(:)-V4(:))),max(abs(Policy1(:)-Policy4(:))),max(abs(StationaryDist1(:)-StationaryDist4(:))))

end

%%
output=struct();

end
