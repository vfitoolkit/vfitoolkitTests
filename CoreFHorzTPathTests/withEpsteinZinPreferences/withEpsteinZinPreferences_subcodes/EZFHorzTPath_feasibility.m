function output=EZFHorzTPath_feasibility(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames)
% The two feasibility conventions of the Epstein-Zin raws (EZ_feasibility_conventions_proposal.md at the
% toolkit root), on the MOVING transition path. The transition-path raws must carry the same four lines
% as the stationary ones; the do-nothing checks of the shape subcodes never reach either convention,
% because their models have strictly signed returns and (away from zero assets at the kappa_j=0 ages) no
% dead states. The trigger models are those of EZFHorz_CrossTests_feasibility in CoreFHorzEZTests.
%
% A: STATES WITH NO FEASIBLE CHOICE (V=-Inf). a_grid extended below the natural borrowing limit with
%    integer-spaced nodes from -12 to -1. The vNM path of the same budget gives the reference dead set.
%    Under negative utility (ezc3=ezc4=-1) a -Inf continuation left untouched by the transform would come
%    out of ezc3*beta*(.) as +Inf, so a choice leading to a dead state would be valued +Inf and chosen.
%    Checks, on the moving path, both layouts: no +Inf and no NaN in V; the -Inf set equals the vNM dead
%    set; the dead set is the same in all three cases; the lowmemory, DC, GI and DC+GI legs give the same
%    dead set and (off GI) the same V on the finite set, and DC+GI the same V as GI.
% B: A RETURN OF EXACTLY ZERO AT A FEASIBLE CHOICE. Subsistence-floor return fns, floorval=0 against
%    floorval=+-1e-10, on the moving path: the same -Inf set, and V equal on the common finite set to the
%    scale of 1e-10, except in consumption units, where the gap is linear in floorval but far above it, so
%    floorval=0 is checked to be the limit instead (a solve at 1e-11 must shrink the gap tenfold). Anti-vacuity: in period 1 of the path the share of states whose chosen consumption is
%    on the floor must be well above zero.

DF=DiscountFactorParamNames;

ReturnFn_cons=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) EZReturnFn_cons_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,agej,Jr,pension);
ReturnFn_posU=@(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension) EZReturnFn_positiveUtils_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_negU=@(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension) EZReturnFn_negativeUtils_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_cons_floor=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension,xbar,floorval) EZReturnFn_cons_nod_z_floor(aprime,a,z,r,w,kappa_j,agej,Jr,pension,xbar,floorval);
ReturnFn_posU_floor=@(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension,xbar,floorval) EZReturnFn_positiveUtils_nod_z_floor(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension,xbar,floorval);
ReturnFn_negU_floor=@(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension,kappa,floorval) EZReturnFn_negativeUtils_nod_z_floor(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension,kappa,floorval);

transpathoptionsL=cell(1,2);
transpathoptionsL{1}.fastOLG=0;
transpathoptionsL{2}.fastOLG=1;
layoutname={'slowOLG','fastOLG'};
casename={'cons units','positive utils','negative utils'};

%% A: states with no feasible choice
n_neg=12;
a_gridNBL=[(-n_neg:1:-1)'; a_grid];
n_aNBL=length(a_gridNBL);
fprintf('A: a_grid has %i integer-spaced nodes from %i to -1 below the usual cubic grid (n_a=%i) \n',n_neg,-n_neg,n_aNBL)
legnames={'lowmemory=1','divide-and-conquer','grid interpolation','DC+GI'};
Vdead=cell(3,2);
for cc=1:3
    vfEZ=struct();
    vfEZ.exoticpreferences='EpsteinZin';
    if cc==1
        vfEZ.EZutils=0; vfEZ.EZriskaversion='ezgamma'; vfEZ.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
        Paramsref=Params; Paramsref.ezsigma=Params.ezgamma; ReturnFnref=ReturnFn_negU;
    elseif cc==2
        vfEZ.EZutils=1; vfEZ.EZpositiveutility=1; vfEZ.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
        Paramsref=Params; ReturnFnref=ReturnFn_posU;
    else
        vfEZ.EZutils=1; vfEZ.EZpositiveutility=0; vfEZ.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
        Paramsref=Params; ReturnFnref=ReturnFn_negU;
    end
    % the vNM reference: the same budget, so the same dead set
    vfN=struct();
    [Vn_s,Policyn_s]=ValueFnIter_Case1_FHorz(0,n_aNBL,n_z,N_j,[],a_gridNBL,z_grid,pi_z,ReturnFnref,Paramsref,DF,[],vfN);
    [V_s,Policy_s]=ValueFnIter_Case1_FHorz(0,n_aNBL,n_z,N_j,[],a_gridNBL,z_grid,pi_z,ReturnFn,Params,DF,[],vfEZ);
    vfGI=vfEZ; vfGI.gridinterplayer=1; vfGI.ngridinterp=5;
    [V_sGI,Policy_sGI]=ValueFnIter_Case1_FHorz(0,n_aNBL,n_z,N_j,[],a_gridNBL,z_grid,pi_z,ReturnFn,Params,DF,[],vfGI);
    for ll=1:2
        VPathn=gather(ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Vn_s, Policyn_s, Paramsref, 0, n_aNBL, n_z, N_j, [], a_gridNBL,z_grid, pi_z, DF, ReturnFnref, transpathoptionsL{ll}, vfN));
        deadN=(VPathn==-Inf);
        fprintf('A [%s, %s] the vNM reference dead set on the path is non-empty, should be WELL ABOVE zero: %i \n',casename{cc},layoutname{ll},sum(deadN(:)))
        VPath1=gather(ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_s, Policy_s, Params, 0, n_aNBL, n_z, N_j, [], a_gridNBL,z_grid, pi_z, DF, ReturnFn, transpathoptionsL{ll}, vfEZ));
        fin=isfinite(VPath1);
        fprintf('A [%s, %s] census of V on the path: finite %i, -Inf %i, +Inf %i, NaN %i of %i \n',casename{cc},layoutname{ll},sum(fin(:)),sum(VPath1(:)==-Inf),sum(VPath1(:)==Inf),sum(isnan(VPath1(:))),numel(VPath1))
        fprintf('A [%s, %s] the finite set is non-empty, this should be one: %d \n',casename{cc},layoutname{ll},sum(fin(:))>0)
        fprintf('A [%s, %s] V has no +Inf, this should be one: %d \n',casename{cc},layoutname{ll},sum(VPath1(:)==Inf)==0)
        fprintf('A [%s, %s] V has no NaN, this should be one: %d \n',casename{cc},layoutname{ll},sum(isnan(VPath1(:)))==0)
        fprintf('A [%s, %s] the -Inf set equals the vNM dead set, this should be one: %d \n',casename{cc},layoutname{ll},isequal(VPath1==-Inf,deadN))
        Vdead{cc,ll}=(VPath1==-Inf);
        VPathGI=[];
        for lg=1:4
            vfL=vfEZ;
            VfinalL=V_s; PolicyfinalL=Policy_s;
            if lg==1
                vfL.lowmemory=1;
            elseif lg==2
                vfL.divideandconquer=1;
            elseif lg==3
                vfL.gridinterplayer=1; vfL.ngridinterp=5;
                VfinalL=V_sGI; PolicyfinalL=Policy_sGI;
            else
                vfL.divideandconquer=1; vfL.gridinterplayer=1; vfL.ngridinterp=5;
                VfinalL=V_sGI; PolicyfinalL=Policy_sGI;
            end
            VPathL=gather(ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, VfinalL, PolicyfinalL, Params, 0, n_aNBL, n_z, N_j, [], a_gridNBL,z_grid, pi_z, DF, ReturnFn, transpathoptionsL{ll}, vfL));
            fprintf('A [%s, %s] %s leg: V has no +Inf and no NaN, this should be one: %d \n',casename{cc},layoutname{ll},legnames{lg},(sum(VPathL(:)==Inf)+sum(isnan(VPathL(:))))==0)
            fprintf('A [%s, %s] %s leg: the -Inf set equals the plain leg''s, this should be one: %d \n',casename{cc},layoutname{ll},legnames{lg},isequal(VPathL==-Inf,VPath1==-Inf))
            if lg<=2
                fprintf('A [%s, %s] %s leg: V on the finite set against the plain leg, this should be zero: %.3e \n',casename{cc},layoutname{ll},legnames{lg},max(abs(VPathL(fin)-VPath1(fin))))
            elseif lg==3
                VPathGI=VPathL;
            else
                finGI=isfinite(VPathGI);
                fprintf('A [%s, %s] DC+GI against GI, V on the finite set, this should be zero: %.3e \n',casename{cc},layoutname{ll},max(abs(VPathL(finGI)-VPathGI(finGI))))
            end
        end
        clear VPathn VPath1 VPathL VPathGI deadN fin
    end
end
for ll=1:2
    fprintf('A [%s] the dead set is the same in all three cases (same budget), this should be one: %d \n',layoutname{ll},isequal(Vdead{1,ll},Vdead{2,ll}) && isequal(Vdead{1,ll},Vdead{3,ll}))
end
clear Vdead

%% B: a return of exactly zero at a feasible choice
Pf=Params;
Pf.xbar=0.5;  % subsistence consumption: at a=0 with kappa_j=0.5 the whole wage is below it, so the floor is FORCED there
Pf.kappa=0.5; % negative utils: u(c)+kappa reaches zero at c=1/kappa=2, so the satiation floor is hit at high assets
epsval=10^(-10);
fprintf('B: xbar=%g (cons units and positive utils), kappa=%g (negative utils), floorval 0 against %g \n',Pf.xbar,Pf.kappa,epsval)
N_z=prod(n_z);
for cc=1:3
    vfEZ=struct();
    vfEZ.exoticpreferences='EpsteinZin';
    if cc==1
        vfEZ.EZutils=0; vfEZ.EZriskaversion='ezgamma'; vfEZ.EZeis='ezphi';
        ReturnFn=ReturnFn_cons_floor;
        epssigned=epsval;
    elseif cc==2
        vfEZ.EZutils=1; vfEZ.EZpositiveutility=1; vfEZ.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU_floor;
        epssigned=epsval;
    else
        vfEZ.EZutils=1; vfEZ.EZpositiveutility=0; vfEZ.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU_floor;
        epssigned=-epsval;
    end
    P0=Pf; P0.floorval=0;
    Pe=Pf; Pe.floorval=epssigned;
    [V0_s,Policy0_s]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn,P0,DF,[],vfEZ);
    [Ve_s,Policye_s]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn,Pe,DF,[],vfEZ);
    for ll=1:2
        VPath0=gather(ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V0_s, Policy0_s, P0, 0, n_a, n_z, N_j, [], a_grid,z_grid, pi_z, DF, ReturnFn, transpathoptionsL{ll}, vfEZ));
        [VPathe,PolicyPathe]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Ve_s, Policye_s, Pe, 0, n_a, n_z, N_j, [], a_grid,z_grid, pi_z, DF, ReturnFn, transpathoptionsL{ll}, vfEZ);
        VPathe=gather(VPathe);
        % anti-vacuity: in period 1 of the path, the floor must be chosen somewhere
        PolicyVals1=gather(PolicyInd2Val_FHorz(PolicyPathe(:,:,:,:,1),0,n_a,n_z,N_j,[],a_grid,vfEZ));
        aprimevals=reshape(PolicyVals1,[n_a,N_z,N_j]);
        cvals=zeros(n_a,N_z,N_j);
        for jj=1:N_j
            if Params.agej(jj)<Params.Jr
                cvals(:,:,jj)=(1+PricePath.r(1))*a_grid+ParamPath.w(1)*Params.kappa_j(jj)*z_grid'-aprimevals(:,:,jj);
            else
                cvals(:,:,jj)=(1+PricePath.r(1))*a_grid+Params.pension-aprimevals(:,:,jj);
            end
        end
        if cc<=2
            floorhit=(cvals>0)&(cvals<Pf.xbar);
        else
            floorhit=(cvals>0)&((cvals.^(1-Params.ezsigma))/(1-Params.ezsigma)+Pf.kappa>=0);
        end
        fprintf('B [%s, %s] share of (a,z,j) in period 1 whose chosen consumption is on the floor, should be WELL ABOVE zero: %2.4f \n',casename{cc},layoutname{ll},mean(floorhit(:)))
        fprintf('B [%s, %s] floorval=0 census of V on the path: finite %i, -Inf %i, +Inf %i, NaN %i, exactly 0: %i of %i \n',casename{cc},layoutname{ll},sum(isfinite(VPath0(:))),sum(VPath0(:)==-Inf),sum(VPath0(:)==Inf),sum(isnan(VPath0(:))),sum(VPath0(:)==0),numel(VPath0))
        fprintf('B [%s, %s] floorval=0: V has no +Inf and no NaN, this should be one: %d \n',casename{cc},layoutname{ll},(sum(VPath0(:)==Inf)+sum(isnan(VPath0(:))))==0)
        fprintf('B [%s, %s] the -Inf set with floorval=0 equals that with floorval=%g, this should be one: %d \n',casename{cc},layoutname{ll},epssigned,isequal(VPath0==-Inf,VPathe==-Inf))
        fin=isfinite(VPath0)&isfinite(VPathe);
        fprintf('B [%s, %s] the common finite set is non-empty, this should be one: %d \n',casename{cc},layoutname{ll},sum(fin(:))>0)
        if cc==1
            % Consumption units (EIS 0.5, risk aversion 3, so negative powers in the aggregator): with floorval=0 any
            % state that can reach a floor state next period has CE=0 and so V exactly 0, while with floorval=eps its V
            % is about eps/(beta*sqrt(p)), p the probability of reaching the floor, which with the tail probabilities of
            % z is far above eps. So the gap is not of the scale of eps (it was 1.4e-6 at eps=1e-10, 2026-11-02: this
            % model hits the floor in 16% of states, the stationary bank's in 0.4%). What must hold is that V(floorval=0)
            % is the eps->0 limit: the aggregator is homogeneous of degree one, so the gap is linear in eps, and a
            % third solve at eps/10 must shrink it exactly tenfold.
            Pe10=Pf; Pe10.floorval=epssigned/10;
            [Ve10_s,Policye10_s]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn,Pe10,DF,[],vfEZ);
            VPathe10=gather(ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, Ve10_s, Policye10_s, Pe10, 0, n_a, n_z, N_j, [], a_grid,z_grid, pi_z, DF, ReturnFn, transpathoptionsL{ll}, vfEZ));
            fin10=fin&isfinite(VPathe10);
            gap1=max(abs(VPath0(fin10)-VPathe(fin10)));
            gap10=max(abs(VPath0(fin10)-VPathe10(fin10)));
            fprintf('B [%s, %s] V with floorval=0 against floorval=%g on the common finite set, DIAGNOSTIC (not of the scale of floorval, see the code): %.3e \n',casename{cc},layoutname{ll},epssigned,gap1)
            fprintf('B [%s, %s] the gap at floorval=%g is not zero (else the tenfold check below is 0/0), should NOT be zero: %.3e \n',casename{cc},layoutname{ll},epssigned/10,gap10)
            fprintf('B [%s, %s] floorval=0 is the limit: the gap shrinks tenfold from floorval=%g to %g, ratio/10-1, this should be close to zero: %.3e \n',casename{cc},layoutname{ll},epssigned,epssigned/10,abs(gap1/(10*gap10)-1))
            clear VPathe10 Ve10_s Policye10_s fin10
        else
            fprintf('B [%s, %s] V with floorval=0 against floorval=%g on the common finite set, this should be zero: %.3e \n',casename{cc},layoutname{ll},epssigned,max(abs(VPath0(fin)-VPathe(fin))))
        end
        clear VPath0 VPathe PolicyPathe PolicyVals1 aprimevals cvals floorhit fin
    end
end

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
