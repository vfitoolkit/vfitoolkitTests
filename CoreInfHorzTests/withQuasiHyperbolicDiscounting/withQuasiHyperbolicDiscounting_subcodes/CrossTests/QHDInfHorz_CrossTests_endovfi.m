function output=QHDInfHorz_CrossTests_endovfi(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test: endogenous VFI (vfoptions.relativevfi) for NAIVE quasi-hyperbolic.
%
% Naive QH first solves the exponential-discounting problem (Valt, Policyalt, discount factor beta), and then takes one
% quasi-hyperbolic step from Valt (discount factor beta0*beta). The first solve is a standard VFI, so endogenous VFI
% applies to it, and is the default from beta>=0.97. It must give the same Valt (in levels, as the QH step reads it) and
% Policyalt, and hence the same V and Policy, as plain VFI.
%
% Paths covered (each is a separate call site for the exponential solve):
%   (a) nod, pure discretization          -> ValueFnIter_InfHorz_[EndogenousVFI_]nod_raw
%   (b) d,   refinement (the default)     -> ValueFnIter_InfHorz_[EndogenousVFI_]nod_raw, inside QuasiHyperbolic_Refine
%   (c) d,   refinement, howardssparse=1  -> ..._[EndogenousVFI_]sparse_nod_raw
%   (d) d,   refinement, howardsgreedy=1  -> ..._[EndogenousVFI_]HowardGreedy_nod_raw
%   (e) d,   pure discretization          -> ValueFnIter_InfHorz_[EndogenousVFI_]raw
%   (f) nod, grid interpolation           -> ValueFnIter_InfHorz_GridInterpLayer (which switches to the GI twins itself)
%   (g) d,   grid interpolation           -> the same, with refinement
% Run at beta=0.99, with tolerance=1e-10: at beta=0.99 each solve is within about beta/(1-beta)*tolerance of the true V,
% which at the default 1e-9 is about 1e-7 and would put the difference of two solves at the identity bar.
% Then: the default (relativevfi not set) at beta=0.99 is identical to the explicit relativevfi=1; and Sophisticated
% still refuses relativevfi=1 (its value is not a standard VFI).

fprintf('\n================ CROSS-TEST: Naive quasi-hyperbolic, endogenous VFI vs plain VFI ================\n');

ReturnFn_nod=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
ReturnFn_d=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
DF=DiscountFactorParamNames;

ParamsE=Params;
ParamsE.beta=0.99;
ParamsE.r=0.8*(1/ParamsE.beta-1); % so that beta*(1+r)<1

vfobase=struct();
vfobase.exoticpreferences='QuasiHyperbolic';
vfobase.quasi_hyperbolic='Naive';
vfobase.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfobase.tolerance=10^(-10);

pathlabels={'(a) nod, pure discretization','(b) d, refinement','(c) d, refinement, sparse','(d) d, refinement, greedy','(e) d, pure discretization','(f) nod, grid interpolation','(g) d, grid interpolation'};
for pp=1:length(pathlabels)
    vfo=vfobase;
    if pp==1 || pp==6
        nd=0; dg=[]; RF=ReturnFn_nod;
    else
        nd=n_d; dg=d_grid; RF=ReturnFn_d;
    end
    if pp==3
        vfo.howardssparse=1;
    elseif pp==4
        vfo.howardsgreedy=1;
    elseif pp==5
        vfo.solnmethod='purediscretization';
    elseif pp>=6
        vfo.gridinterplayer=1; vfo.ngridinterp=5;
    end
    vfo0=vfo; vfo0.relativevfi=0;
    [V0,Policy0,Valt0,Policyalt0]=ValueFnIter_InfHorz(nd,n_a,n_z,dg,a_grid,z_grid,pi_z,RF,ParamsE,DF,[],vfo0);
    vfo1=vfo; vfo1.relativevfi=1;
    [V1,Policy1,Valt1,Policyalt1]=ValueFnIter_InfHorz(nd,n_a,n_z,dg,a_grid,z_grid,pi_z,RF,ParamsE,DF,[],vfo1);
    fprintf('QH Naive %-32s endoVFI vs plain, V         should be ~0: %.3e \n',pathlabels{pp},max(abs(V0(:)-V1(:))));
    fprintf('QH Naive %-32s endoVFI vs plain, Policy    should be  0: %.3e \n',pathlabels{pp},max(abs(Policy0(:)-Policy1(:))));
    fprintf('QH Naive %-32s endoVFI vs plain, Valt      should be ~0: %.3e \n',pathlabels{pp},max(abs(Valt0(:)-Valt1(:))));
    fprintf('QH Naive %-32s endoVFI vs plain, Policyalt should be  0: %.3e \n',pathlabels{pp},max(abs(Policyalt0(:)-Policyalt1(:))));
end

% The default at beta=0.99 is endogenous VFI with 40 Howards steps (indexed Howards, as at these grid sizes without GI)
for pp=1:2
    if pp==1
        nd=0; dg=[]; RF=ReturnFn_nod;
    else
        nd=n_d; dg=d_grid; RF=ReturnFn_d;
    end
    vfo=vfobase;
    [VD,PolicyD,ValtD]=ValueFnIter_InfHorz(nd,n_a,n_z,dg,a_grid,z_grid,pi_z,RF,ParamsE,DF,[],vfo);
    vfo=vfobase; vfo.relativevfi=1; vfo.howards=40; vfo.howardssparse=0;
    [VX,PolicyX,ValtX]=ValueFnIter_InfHorz(nd,n_a,n_z,dg,a_grid,z_grid,pi_z,RF,ParamsE,DF,[],vfo);
    fprintf('QH Naive %-32s default at beta=0.99 vs explicit endoVFI, V      should be  0: %.3e \n',pathlabels{pp},max(abs(VD(:)-VX(:))));
    fprintf('QH Naive %-32s default at beta=0.99 vs explicit endoVFI, Policy should be  0: %.3e \n',pathlabels{pp},max(abs(PolicyD(:)-PolicyX(:))));
    fprintf('QH Naive %-32s default at beta=0.99 vs explicit endoVFI, Valt   should be  0: %.3e \n',pathlabels{pp},max(abs(ValtD(:)-ValtX(:))));
end

% Sophisticated still refuses an explicit relativevfi=1
vfo=vfobase; vfo.quasi_hyperbolic='Sophisticated'; vfo.relativevfi=1;
ok=0; msg='';
try
    [~,~]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,ParamsE,DF,[],vfo);
catch ME
    msg=ME.message;
    ok=double(contains(msg,'relativevfi=1 cannot be used with vfoptions.exoticpreferences'));
end
fprintf('QH Sophisticated with relativevfi=1 is refused,                      should be one: %i  (%s)\n',ok,msg);

output=struct();

end
