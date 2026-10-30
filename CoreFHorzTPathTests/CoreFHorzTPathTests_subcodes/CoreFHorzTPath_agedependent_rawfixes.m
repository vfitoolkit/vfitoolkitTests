function output=CoreFHorzTPath_agedependent_rawfixes(T,PricePath,ParamPath,n_d,n_a,n_z,n_e,N_j,d_grid,a_grid,z_grid,pi_z,e_grid,pi_e,Params,DiscountFactorParamNames,vfoptionsbaseline)
% Checks for four slowOLG/fastOLG transition-path raw bugs fixed 2026-10-28. No figures, no general eqm.
% Every other part of this bank uses e and z grids and pi_e that do NOT vary with age, and under those
% the first two bugs give the right numbers by accident: that is how they survived a 956-check run. So
% here the shocks VARY BY AGE, and each section first checks that the variation actually matters
% (a 'should NOT be zero'), so that the identity after it cannot be vacuous.
%   (A) the noz_e GI raws took the e expectation with pi_e_J(:,jj) instead of jj+1 (Reading A). Oracle:
%       the DC+GI raw, which was always right, so GI against DC+GI must agree.
%   (B) GI1_nod_e_raw under lowmemory=1/2 read e_gridvals_J and z_gridvals_J without the age index.
%       lowmemory=1 and =2 must agree with lowmemory=0.
%   (C) the with-d DC1 raws errored with a single d point: the DC index arrays collapse to vectors at
%       N_d=1, and a_grid(vector) is a column whatever the index's orientation. DC must agree with no DC
%       at N_d=1, for every shock case and both layouts.
%   (D) the fastOLG non-DC noz, e and noz_e raws ignore vfoptions.EVpre; EVpre=1 is now refused there.

transpathoptionsslow.fastOLG=0;
Params.w=1;
Params.r=0.05;

% Age-dependent shocks. The e grid scales with age, and pi_e moves from pi_e towards uniform with age.
e_grid_J=e_grid.*linspace(0.8,1.2,N_j);                            % [n_e,N_j]
wJ=linspace(0,0.8,N_j);
pi_e_J=pi_e.*(1-wJ)+(1/n_e)*wJ;                                     % [n_e,N_j], each column sums to one
pi_e_Jrev=pi_e_J(:,end:-1:1);                                       % the same columns in reverse age order
z_grid_J=z_grid.*linspace(0.9,1.1,N_j);                             % [n_z,N_j]

%% (A) noz_e GI: pi_e timing
for dd=0:1
    if dd==0
        n_d_cc=0; d_grid_cc=[];
        ReturnFn=@(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension);
        Policy_final=ones([3,n_a,n_e,N_j],'gpuArray');
        lab='age-dependent e (A, nod noz_e)';
    else
        n_d_cc=n_d; d_grid_cc=d_grid;
        ReturnFn=@(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
        Policy_final=ones([4,n_a,n_e,N_j],'gpuArray');
        lab='age-dependent e (A, d noz_e)';
    end
    V_final=zeros([n_a,n_e,N_j],'gpuArray');
    vfoptions=struct();
    vfoptions.n_e=n_e; vfoptions.e_grid=e_grid_J; vfoptions.pi_e=pi_e_J;
    vfoptions.gridinterplayer=1; vfoptions.ngridinterp=5;
    VPathGI=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d_cc, n_a, 0, N_j, d_grid_cc, a_grid,[], [], DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions);
    vfoptionsR=vfoptions; vfoptionsR.pi_e=pi_e_Jrev;
    VPathGIrev=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d_cc, n_a, 0, N_j, d_grid_cc, a_grid,[], [], DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptionsR);
    fprintf('%s: the age profile of pi_e matters (reversing it moves V), should NOT be zero: %.3e \n',lab,max(abs(VPathGI(:)-VPathGIrev(:))))
    vfoptionsD=vfoptions; vfoptionsD.divideandconquer=1;
    VPathDCGI=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d_cc, n_a, 0, N_j, d_grid_cc, a_grid,[], [], DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptionsD);
    fprintf('%s: GI against DC+GI, V, this should be zero: %.3e \n',lab,max(abs(VPathGI(:)-VPathDCGI(:))))
    clear VPathGI VPathGIrev VPathDCGI
end

%% (B) GI1_nod_e lowmemory: the age index on the shock grids
ReturnFn=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
V_final=zeros([n_a,n_z,n_e,N_j],'gpuArray');
Policy_final=ones([3,n_a,n_z,n_e,N_j],'gpuArray');
vfoptions=struct();
vfoptions.n_e=n_e; vfoptions.e_grid=e_grid_J; vfoptions.pi_e=pi_e_J;
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=5;
lab='age-dependent e and z (B, nod z e, GI)';
VPath0=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, 0, n_a, n_z, N_j, [], a_grid,z_grid_J, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptions);
vfoptionsC=vfoptions; vfoptionsC.e_grid=e_grid.*ones(1,N_j);
VPathC=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, 0, n_a, n_z, N_j, [], a_grid,z_grid_J, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptionsC);
fprintf('%s: the age profile of the e grid matters, should NOT be zero: %.3e \n',lab,max(abs(VPath0(:)-VPathC(:))))
for lm=1:2
    vfoptionsM=vfoptions; vfoptionsM.lowmemory=lm;
    VPathM=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, 0, n_a, n_z, N_j, [], a_grid,z_grid_J, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsslow, vfoptionsM);
    fprintf('%s: lowmemory=%i against lowmemory=0, V, this should be zero: %.3e \n',lab,lm,max(abs(VPathM(:)-VPath0(:))))
end
clear VPath0 VPathC VPathM

%% (C) Divide-and-conquer with a single d point
% With N_d=1 the DC index arrays (loweredge, aprimeindexes) can collapse to vectors, and indexing the
% column a_grid by a vector returns a COLUMN whatever the index's orientation, so the window grid handed to
% the return-matrix builder had the wrong shape ('Number of elements must not change'). Every shock case
% and both layouts, lowmemory 0 and 1 (the lowmemory branches index per z/e, where it collapses most).
d_grid1=d_grid(ceil(length(d_grid)/2));
layoutname={'slowOLG','fastOLG'};
shockname={'z','noz','e','noz_e'};
for ll=1:2
    transpathoptionsC.fastOLG=ll-1;
    for ss=1:4
        vfoptions=struct();
        if ss==1
            ReturnFn=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
            V_final=zeros([n_a,n_z,N_j],'gpuArray'); Policy_final=ones([2,n_a,n_z,N_j],'gpuArray');
            n_z_cc=n_z; z_grid_cc=z_grid; pi_z_cc=pi_z;
        elseif ss==2
            ReturnFn=@(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_noe_nosemiz(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
            V_final=zeros([n_a,N_j],'gpuArray'); Policy_final=ones([2,n_a,N_j],'gpuArray');
            n_z_cc=0; z_grid_cc=[]; pi_z_cc=[];
        elseif ss==3
            ReturnFn=@(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
            V_final=zeros([n_a,n_z,n_e,N_j],'gpuArray'); Policy_final=ones([2,n_a,n_z,n_e,N_j],'gpuArray');
            n_z_cc=n_z; z_grid_cc=z_grid; pi_z_cc=pi_z;
            vfoptions.n_e=n_e; vfoptions.e_grid=e_grid; vfoptions.pi_e=pi_e;
        else
            ReturnFn=@(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
            V_final=zeros([n_a,n_e,N_j],'gpuArray'); Policy_final=ones([2,n_a,n_e,N_j],'gpuArray');
            n_z_cc=0; z_grid_cc=[]; pi_z_cc=[];
            vfoptions.n_e=n_e; vfoptions.e_grid=e_grid; vfoptions.pi_e=pi_e;
        end
        lab=['single d point (C, ',layoutname{ll},', ',shockname{ss},', DC)'];
        VPath0=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, 1, n_a, n_z_cc, N_j, d_grid1, a_grid,z_grid_cc, pi_z_cc, DiscountFactorParamNames, ReturnFn, transpathoptionsC, vfoptions);
        for lm=0:1
            vfoptionsD=vfoptions; vfoptionsD.divideandconquer=1; vfoptionsD.lowmemory=lm;
            try
                VPathD=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, 1, n_a, n_z_cc, N_j, d_grid1, a_grid,z_grid_cc, pi_z_cc, DiscountFactorParamNames, ReturnFn, transpathoptionsC, vfoptionsD);
                fprintf('%s: lowmemory=%i, DC against no DC, V, this should be zero: %.3e \n',lab,lm,max(abs(VPathD(:)-VPath0(:))))
            catch ME
                fprintf('%s: lowmemory=%i, DC errored (%s), this should be zero: %.3e \n',lab,lm,ME.message,1)
            end
        end
    end
end
clear VPath0 VPathD

%% (C2) The same, through the grid interpolation layer, quasi-hyperbolic discounting and a semi-exogenous state
% The DC+GI, QH and SemiExo divide-and-conquer raws have the same index arrays and the same lookups, so the
% same fix. Each case: divide-and-conquer against no divide-and-conquer at a single d point (a single d1
% point, with the binary d2, for semiz), both layouts. GI with every shock case and lowmemory 0 and 1; QH
% Naive and Sophisticated with and without GI; semiz with and without z, with and without GI.
for ll=1:2
    transpathoptionsC.fastOLG=ll-1;
    % (C2a) grid interpolation layer
    for ss=1:4
        vfoptions=struct(); vfoptions.gridinterplayer=1; vfoptions.ngridinterp=5;
        if ss==1
            ReturnFn=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
            V_final=zeros([n_a,n_z,N_j],'gpuArray'); Policy_final=ones([4,n_a,n_z,N_j],'gpuArray');
            n_z_cc=n_z; z_grid_cc=z_grid; pi_z_cc=pi_z;
        elseif ss==2
            ReturnFn=@(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_noe_nosemiz(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
            V_final=zeros([n_a,N_j],'gpuArray'); Policy_final=ones([4,n_a,N_j],'gpuArray');
            n_z_cc=0; z_grid_cc=[]; pi_z_cc=[];
        elseif ss==3
            ReturnFn=@(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
            V_final=zeros([n_a,n_z,n_e,N_j],'gpuArray'); Policy_final=ones([4,n_a,n_z,n_e,N_j],'gpuArray');
            n_z_cc=n_z; z_grid_cc=z_grid; pi_z_cc=pi_z;
            vfoptions.n_e=n_e; vfoptions.e_grid=e_grid; vfoptions.pi_e=pi_e;
        else
            ReturnFn=@(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
            V_final=zeros([n_a,n_e,N_j],'gpuArray'); Policy_final=ones([4,n_a,n_e,N_j],'gpuArray');
            n_z_cc=0; z_grid_cc=[]; pi_z_cc=[];
            vfoptions.n_e=n_e; vfoptions.e_grid=e_grid; vfoptions.pi_e=pi_e;
        end
        lab=['single d point (C2a, ',layoutname{ll},', ',shockname{ss},', DC+GI)'];
        VPath0=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, 1, n_a, n_z_cc, N_j, d_grid1, a_grid,z_grid_cc, pi_z_cc, DiscountFactorParamNames, ReturnFn, transpathoptionsC, vfoptions);
        for lm=0:1
            vfoptionsD=vfoptions; vfoptionsD.divideandconquer=1; vfoptionsD.lowmemory=lm;
            try
                VPathD=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, 1, n_a, n_z_cc, N_j, d_grid1, a_grid,z_grid_cc, pi_z_cc, DiscountFactorParamNames, ReturnFn, transpathoptionsC, vfoptionsD);
                fprintf('%s: lowmemory=%i, DC+GI against GI, V, this should be zero: %.3e \n',lab,lm,max(abs(VPathD(:)-VPath0(:))))
            catch ME
                fprintf('%s: lowmemory=%i, DC+GI errored (%s), this should be zero: %.3e \n',lab,lm,ME.message,1)
            end
        end
    end
    % (C2b) quasi-hyperbolic discounting, Naive and Sophisticated, with z (lowmemory 0 and 1), without z, and with e but no z
    ParamsQH=Params; ParamsQH.beta0=0.9;
    qhname={'Naive','Sophisticated'};
    for qq=1:2
        for ss=[1,2,4]
            for gi=0:1
                vfoptions=struct(); vfoptions.exoticpreferences='QuasiHyperbolic'; vfoptions.quasi_hyperbolic=qhname{qq}; vfoptions.QHadditionaldiscount='beta0';
                if gi==1
                    vfoptions.gridinterplayer=1; vfoptions.ngridinterp=5;
                end
                if ss==1
                    ReturnFn=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
                    V_final=zeros([n_a,n_z,N_j],'gpuArray'); Policy_final=ones([2+2*gi,n_a,n_z,N_j],'gpuArray');
                    n_z_cc=n_z; z_grid_cc=z_grid; pi_z_cc=pi_z;
                elseif ss==2
                    ReturnFn=@(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_noe_nosemiz(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
                    V_final=zeros([n_a,N_j],'gpuArray'); Policy_final=ones([2+2*gi,n_a,N_j],'gpuArray');
                    n_z_cc=0; z_grid_cc=[]; pi_z_cc=[];
                else
                    ReturnFn=@(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
                    V_final=zeros([n_a,n_e,N_j],'gpuArray'); Policy_final=ones([2+2*gi,n_a,n_e,N_j],'gpuArray');
                    n_z_cc=0; z_grid_cc=[]; pi_z_cc=[];
                    vfoptions.n_e=n_e; vfoptions.e_grid=e_grid; vfoptions.pi_e=pi_e;
                end
                VPath0=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, ParamsQH, 1, n_a, n_z_cc, N_j, d_grid1, a_grid,z_grid_cc, pi_z_cc, DiscountFactorParamNames, ReturnFn, transpathoptionsC, vfoptions);
                for lm=0:(ss==1) % lowmemory=1 loops over z, so only with z
                    lab=sprintf('single d point (C2b, %s, QH %s, %s, GI%i, DC, lowmemory=%i)',layoutname{ll},qhname{qq},shockname{ss},gi,lm);
                    vfoptionsD=vfoptions; vfoptionsD.divideandconquer=1; vfoptionsD.lowmemory=lm;
                    try
                        VPathD=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, ParamsQH, 1, n_a, n_z_cc, N_j, d_grid1, a_grid,z_grid_cc, pi_z_cc, DiscountFactorParamNames, ReturnFn, transpathoptionsC, vfoptionsD);
                        fprintf('%s: DC against no DC, V, this should be zero: %.3e \n',lab,max(abs(VPathD(:)-VPath0(:))))
                    catch ME
                        fprintf('%s: DC errored (%s), this should be zero: %.3e \n',lab,ME.message,1)
                    end
                end
            end
        end
    end
    % (C2c) a semi-exogenous state, with a single d1 point and the binary d2 that drives it
    d_grid_semiz1=[d_grid1; 0; 1];
    for ss=1:2
        for gi=0:1
            vfoptions=struct();
            vfoptions.n_semiz=vfoptionsbaseline.n_semiz; vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid; vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
            if gi==1
                vfoptions.gridinterplayer=1; vfoptions.ngridinterp=5;
            end
            if ss==1
                ReturnFn=@(d1,d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_z_noe_semiz(d1,d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost);
                V_final=zeros([n_a,vfoptions.n_semiz,n_z,N_j],'gpuArray'); Policy_final=ones([3+2*gi,n_a,vfoptions.n_semiz,n_z,N_j],'gpuArray');
                n_z_cc=n_z; z_grid_cc=z_grid; pi_z_cc=pi_z; lab0='semiz, z';
            else
                ReturnFn=@(d1,d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_noz_noe_semiz(d1,d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost);
                V_final=zeros([n_a,vfoptions.n_semiz,N_j],'gpuArray'); Policy_final=ones([3+2*gi,n_a,vfoptions.n_semiz,N_j],'gpuArray');
                n_z_cc=0; z_grid_cc=[]; pi_z_cc=[]; lab0='semiz, noz';
            end
            lab=sprintf('single d1 point (C2c, %s, %s, GI%i, DC)',layoutname{ll},lab0,gi);
            VPath0=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, [1,2], n_a, n_z_cc, N_j, d_grid_semiz1, a_grid,z_grid_cc, pi_z_cc, DiscountFactorParamNames, ReturnFn, transpathoptionsC, vfoptions);
            vfoptionsD=vfoptions; vfoptionsD.divideandconquer=1;
            try
                VPathD=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, [1,2], n_a, n_z_cc, N_j, d_grid_semiz1, a_grid,z_grid_cc, pi_z_cc, DiscountFactorParamNames, ReturnFn, transpathoptionsC, vfoptionsD);
                fprintf('%s: DC against no DC, V, this should be zero: %.3e \n',lab,max(abs(VPathD(:)-VPath0(:))))
            catch ME
                fprintf('%s: DC errored (%s), this should be zero: %.3e \n',lab,ME.message,1)
            end
        end
    end
end
clear VPath0 VPathD

%% (D) EVpre=1 is refused by the fastOLG dispatchers whose non-DC raws ignore it
% The refusal is the first thing these dispatchers do, so the other inputs are never read.
vfo=struct(); vfo.EVpre=1; vfo.divideandconquer=0;
err=zeros(1,3);
try
    ValueFnIter_FHorz_TPath_SingleStep_fastOLG_noz([],n_d,n_a,N_j,[],[],[],struct(),{},{},[],vfo);
catch ME
    err(1)=contains(ME.message,'vfoptions.EVpre=1 (Matched Expectations Path) is only implemented with vfoptions.divideandconquer=1');
end
try
    ValueFnIter_FHorz_TPath_SingleStep_fastOLG_e([],n_d,n_a,n_z,n_e,N_j,[],[],[],[],[],[],[],struct(),{},{},[],vfo);
catch ME
    err(2)=contains(ME.message,'vfoptions.EVpre=1 (Matched Expectations Path) is only implemented with vfoptions.divideandconquer=1');
end
try
    ValueFnIter_FHorz_TPath_SingleStep_fastOLG_noz_e([],n_d,n_a,n_e,N_j,[],[],[],[],[],struct(),{},{},[],vfo);
catch ME
    err(3)=contains(ME.message,'vfoptions.EVpre=1 (Matched Expectations Path) is only implemented with vfoptions.divideandconquer=1');
end
fprintf('EVpre (D): fastOLG noz dispatcher refuses EVpre=1 without DC by name, should be one: %i \n',err(1))
fprintf('EVpre (D): fastOLG e dispatcher refuses EVpre=1 without DC by name, should be one: %i \n',err(2))
fprintf('EVpre (D): fastOLG noz_e dispatcher refuses EVpre=1 without DC by name, should be one: %i \n',err(3))

output=1;

end
