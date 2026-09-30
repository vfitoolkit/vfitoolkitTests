function output=CoreFHorzTPathAlgo_nod_z_noe_nosemiz(T,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,figure_c)
% Compares the transition-path algorithms on one FHorz model: GEnewprice=1 (quasi-Newton with Broyden)
% and GEnewprice=2 (Anderson acceleration) against GEnewprice=3 (shooting), the last with and without
% the additionalfactor ramp. All of them solve the same null-reform general eqm path from the same
% bumped starting guess, and all are timed. The FHorz counterpart of CoreInfHorzTPathAlgo_nod_z_noe_nosemiz,
% on the life-cycle model of CoreFHorzTPath_nod_z_noe_nosemiz.
%
% Plus the checks that are particular to FHorz, or to how these algorithms were brought to it:
%   (G) the triangular FullJacobian against the plain loop, with and without fastOLG
%   (H) the two routes for building the new price path (transpathoptions.updatepert=0 and =1)
%   (I) fastOLG=0 against fastOLG=1, for the shooting path and for the fake-news Jacobian
%   (J) the fake-news (LudwigSSJ) Jacobian against the brute-force FullJacobian

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];

% zero assets, mid point for the shock
jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

%% General equilibrium: solve the stationary GE first, then run a null-reform transition path
% A Cobb-Douglas firm supplies BOTH prices from its first order conditions:
%     r = firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha)) - firmdelta   (MPK less depreciation)
%     w = (1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha))                (MPL)
% so w is a general eqm price here, not the ParamPath entry it is in CoreFHorzTPathTests.
% K is aggregate assets and N aggregate efficiency units of labour (kappa_j*z, with the age weights).
% The GeneralEqmEqns are in residual form, which is what all three algorithms want.
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

FnsToEvaluateGE.K=@(aprime,a,z) a;
FnsToEvaluateGE.N=@(aprime,a,z,kappa_j) kappa_j*z;

GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

heteroagentoptionsGE=struct(); % default fminalgo, and no constraints on r or w

[p_eqm,GEcondns]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluateGE,GeneralEqmEqnsGE,Params,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptionsGE,simoptions,vfoptions);
fprintf('Stationary GE: r=%2.8f, w=%2.8f \n',p_eqm.r,p_eqm.w)
fprintf('Stationary GE conditions (these should be close to zero): CapitalMarket=%.3e, LabourMarket=%.3e \n',GEcondns.CapitalMarket,GEcondns.LabourMarket)

% The path needs the EQUILIBRIUM V_final and initial dist, else the null reform is not exact
ParamsGE=Params;
ParamsGE.r=p_eqm.r;
ParamsGE.w=p_eqm.w;
[V_finalGE,Policy_finalGE]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions);
AgentDist_initialGE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_finalGE,n_d,n_a,n_z,N_j,pi_z,ParamsGE,simoptions);

PricePathGE.r=p_eqm.r*ones(1,T);
PricePathGE.w=p_eqm.w*ones(1,T);
ParamPathGE.sigma=Params.sigma*ones(1,T); % constant: nothing actually changes, this is a null reform

transpathoptionsGE=transpathoptionsbaseline; % fastOLG=1
transpathoptionsGE.maxiter=25;
transpathoptionsGE.verbose=0;
% GEnewprice=3: treat the GeneralEqmEqns as residuals and update each price by a damped fraction of its
% own residual. Both residuals are (price - its firm FOC), so a POSITIVE residual means the price is too
% high and must come down: hence add=0 (subtract) for both.
transpathoptionsGE.GEnewprice=3;
transpathoptionsGE.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; ...
                                            'LabourMarket','w',0,0.1};

%% A bumped starting guess, so all the algorithms have to actually travel
Tbump=floor(2*T/3);
bumpr=0.01+0.04*linspace(0,1,Tbump);
bumpw=0.05-0.04*linspace(0,1,Tbump);
PricePathBumped=PricePathGE;
PricePathBumped.r(1:Tbump)=PricePathGE.r(1:Tbump).*(1+bumpr);
PricePathBumped.w(1:Tbump)=PricePathGE.w(1:Tbump).*(1+bumpw);
fprintf('Bumped starting guess, max initial deviation from p_eqm: r %2.8f, w %2.8f \n',max(abs(PricePathBumped.r-p_eqm.r)),max(abs(PricePathBumped.w-p_eqm.w)))

%% (B) GEnewprice=3: the shooting algorithm, no ramp
transpathoptionsB=transpathoptionsGE;
transpathoptionsB.maxiter=250;
tic;
[PricePathB,GEcondnPathB]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsB, simoptions, vfoptions);
time_B=toc;

%% (C) GEnewprice=3 again, this time with the additionalfactor ramp
transpathoptionsC=transpathoptionsB;
transpathoptionsC.GEnewprice3.additionalfactor=[3,10,30];
tic;
[PricePathC,GEcondnPathC]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsC, simoptions, vfoptions);
time_C=toc;

%% (F) GEnewprice=2: Anderson acceleration of the same shooting map as (B)
% howtoupdate is taken from (B) verbatim: it is the same map, all that differs is how the iterates are
% combined.
transpathoptionsF=rmfield(transpathoptionsB,'GEnewprice3');
transpathoptionsF.GEnewprice=2;
transpathoptionsF.GEnewprice2.howtoupdate=transpathoptionsB.GEnewprice3.howtoupdate;

% (F1) memory=0 empties the history as fast as it fills it, so every step is a plain shooting step and
% Anderson is then exactly the shooting algorithm. This checks the GEnewprice=2 plumbing (howtoupdate
% parsed through GEnewprice2, the whole-path update, the flattening of the path, the terminal row
% staying put), independently of whether Anderson accelerates anything. Both run a fixed 10 iterations,
% short of convergence, so each does the same number of evaluate-then-step cycles.
transpathoptionsF1=transpathoptionsF;
transpathoptionsF1.maxiter=10;
transpathoptionsF1.anderson.memory=0;
[PricePathF0,~]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsF1, simoptions, vfoptions);
transpathoptionsB1=transpathoptionsB;
transpathoptionsB1.maxiter=10;
[PricePathB1,~]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsB1, simoptions, vfoptions);

% (F2) And now with a history to work with. The same five configurations as the InfHorz bank.
memopts=[3,5,10,5,5];
safeopts=[0,0,0,1,0];
typeopts={'II','II','II','II','I'};
nanderson=length(memopts);
time_F=zeros(1,nanderson); distF=zeros(1,nanderson);
devF_r=zeros(1,nanderson); devF_w=zeros(1,nanderson);
andersonname=cell(1,nanderson);
for mm=1:nanderson
    andersonname{mm}=sprintf('memory=%2i, safeguard=%i, type=%s',memopts(mm),safeopts(mm),typeopts{mm});
    transpathoptionsF.anderson.memory=memopts(mm);
    transpathoptionsF.anderson.safeguard=safeopts(mm);
    transpathoptionsF.anderson.type=typeopts{mm};
    fprintf('--- Anderson acceleration, %s --- \n',andersonname{mm})
    tic;
    [PricePathF,GEcondnPathF]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsF, simoptions, vfoptions);
    time_F(mm)=toc;
    distF(mm)=max(sqrt(GEcondnPathF.CapitalMarket.^2+GEcondnPathF.LabourMarket.^2));
    devF_r(mm)=max(abs(PricePathF.r-PricePathB.r));
    devF_w(mm)=max(abs(PricePathF.w-PricePathB.w));
end

%% (A) GEnewprice=1: quasi-Newton on the whole price path, for every Jacobian method and every way of regularising the step
% Run last of the solves, after the cheap shooting ones, so their answers are already in the diary if a
% Newton solve fails.
jacopts={'LudwigPath','LudwigStationary','FullJacobian','LudwigSSJ'};
regopts={'minimumnorm','TikhonovRegularisation','stepcap'};
ncombo=length(jacopts)*length(regopts);
distA=zeros(1,ncombo); time_A=zeros(1,ncombo); devA_r=zeros(1,ncombo); devA_w=zeros(1,ncombo); comboname=cell(1,ncombo);
cc=0;
for jj=1:length(jacopts)
    for rr=1:length(regopts)
        cc=cc+1;
        comboname{cc}=[jacopts{jj},' / ',regopts{rr}];
        transpathoptionsA=transpathoptionsGE;
        transpathoptionsA.maxiter=100;
        transpathoptionsA.GEnewprice=1;
        transpathoptionsA.verbose=1; % quasi-Newton prints a few lines per iteration. Do NOT turn this on
        % for the shooting solves: with GEnewprice=3 verbose dumps the whole price path every iteration.
        transpathoptionsA.GEnewprice1.Jacobianmethod=jacopts{jj};
        transpathoptionsA.GEnewprice1.BroydenRegularisation=regopts{rr};
        fprintf('--- quasi-Newton, Jacobianmethod=%s, BroydenRegularisation=%s --- \n',jacopts{jj},regopts{rr})
        tic;
        [PricePathA,GEcondnPathA]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsA, simoptions, vfoptions);
        time_A(cc)=toc;
        distA(cc)=max(sqrt(GEcondnPathA.CapitalMarket.^2+GEcondnPathA.LabourMarket.^2));
        devA_r(cc)=max(abs(PricePathA.r-PricePathB.r));
        devA_w(cc)=max(abs(PricePathA.w-PricePathB.w));
    end
end

%% Each algorithm must reach its own stated tolerance
distB=max(sqrt(GEcondnPathB.CapitalMarket.^2+GEcondnPathB.LabourMarket.^2));
distC=max(sqrt(GEcondnPathC.CapitalMarket.^2+GEcondnPathC.LabourMarket.^2));
for cc=1:ncombo
    fprintf('quasi-Newton (%-34s), GE condn distance: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',comboname{cc},distA(cc),1e-4)
end
fprintf('Shooting,                              GE condn distance: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',distB,1e-4)
fprintf('Shooting+ramp,                         GE condn distance: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',distC,1e-4)
for mm=1:nanderson
    fprintf('Anderson (%-31s), GE condn distance: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',andersonname{mm},distF(mm),1e-4)
end

%% And they must agree on where the equilibrium path is. Not exactly: each stops the first iteration its
% own residual drops under the tolerance, so they land at different points inside it.
for cc=1:ncombo
    fprintf('quasi-Newton (%-34s) vs shooting, same eqm path, should be small, r: %2.10f, w: %2.10f \n',comboname{cc},devA_r(cc),devA_w(cc))
end
fprintf('Shooting with/without ramp, same eqm path, should be small, r: %2.10f, w: %2.10f \n',max(abs(PricePathC.r-PricePathB.r)),max(abs(PricePathC.w-PricePathB.w)))
fprintf('Anderson with memory=0 vs shooting, both at maxiter=10, the same algorithm, should be zero, r: %.3e, w: %.3e \n',max(abs(PricePathF0.r-PricePathB1.r)),max(abs(PricePathF0.w-PricePathB1.w)))
for mm=1:nanderson
    fprintf('Anderson (%-31s) vs shooting, same eqm path, should be small, r: %2.10f, w: %2.10f \n',andersonname{mm},devF_r(mm),devF_w(mm))
end
for cc=1:ncombo
    fprintf('Runtime, quasi-Newton (%-34s): %2.2f seconds \n',comboname{cc},time_A(cc))
end
fprintf('Runtime, shooting %2.2f seconds, shooting+ramp %2.2f seconds \n',time_B,time_C)
for mm=1:nanderson
    fprintf('Runtime, Anderson (%-31s): %2.2f seconds \n',andersonname{mm},time_F(mm))
end

%% (D) The order of the howtoupdate rows must not matter
% setupGEnewprice3_shooting reorders the rows into PricePath order itself, so giving them in the other
% order has to produce the same path. One iteration is not enough on this TransitionPath: shooting
% returns the path it evaluated, which after one iteration is the starting guess whatever the rows say,
% so two iterations are run and the second path is the first update.
transpathoptionsD=transpathoptionsGE;
transpathoptionsD.maxiter=2;
PricePathD1=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsD, simoptions, vfoptions);
transpathoptionsD.GEnewprice3.howtoupdate={'LabourMarket','w',0,0.1; ...
                                          'CapitalMarket','r',0,0.1}; % the other order
PricePathD2=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsD, simoptions, vfoptions);
fprintf('howtoupdate rows in price order vs eqn order, two iterations, this should be zero, r: %.3e, w: %.3e \n',max(abs(PricePathD2.r-PricePathD1.r)),max(abs(PricePathD2.w-PricePathD1.w)))

%% (H) The two routes for building the new price path must agree
% updatepert=0 (the default) evaluates the general eqm conditions for every period and then updates the
% whole path at once (updatePricePathNew_TPath_T); updatepert=1 updates each period inside the loop over
% t (updatePricePathNew_TPath_tt, the route this TransitionPath always used before). The shooting
% algorithm gives the same answer either way, to the last bit, since each price is updated from its own
% period's condition.
transpathoptionsH=transpathoptionsB1; % 10 iterations of shooting, short of convergence
transpathoptionsH.updatepert=0;
PricePathH0=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsH, simoptions, vfoptions);
transpathoptionsH.updatepert=1;
PricePathH1=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsH, simoptions, vfoptions);
fprintf('updatepert=0 vs updatepert=1, shooting at maxiter=10, this should be zero, r: %.3e, w: %.3e \n',max(abs(PricePathH0.r-PricePathH1.r)),max(abs(PricePathH0.w-PricePathH1.w)))

%% (I) fastOLG=0 vs fastOLG=1
% The single path solve (TransitionPath_FHorz_singlepathiter) has a branch for each, and so does the
% fake-news Jacobian's backward pass. The value fn and agent dist commands agree exactly between the two
% (CoreFHorzTPathTests checks it), so these must too.
transpathoptionsI=transpathoptionsB1; % 10 iterations of shooting
transpathoptionsI.fastOLG=0;
PricePathI0=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsI, simoptions, vfoptions);
fprintf('Shooting at maxiter=10, fastOLG=0 vs fastOLG=1, this should be zero, r: %.3e, w: %.3e \n',max(abs(PricePathI0.r-PricePathB1.r)),max(abs(PricePathI0.w-PricePathB1.w)))
% One quasi-Newton step with the fake-news Jacobian: one Jacobian built, one step taken, so any
% difference in the Jacobian shows up in the price path
transpathoptionsI=transpathoptionsGE;
transpathoptionsI.GEnewprice=1;
transpathoptionsI.maxiter=1;
transpathoptionsI.GEnewprice1.Jacobianmethod='LudwigSSJ';
transpathoptionsI.fastOLG=1;
PricePathISSJ1=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsI, simoptions, vfoptions);
transpathoptionsI.fastOLG=0;
PricePathISSJ0=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsI, simoptions, vfoptions);
fprintf('LudwigSSJ one quasi-Newton step, fastOLG=0 vs fastOLG=1, this should be zero, r: %.3e, w: %.3e \n',max(abs(PricePathISSJ0.r-PricePathISSJ1.r)),max(abs(PricePathISSJ0.w-PricePathISSJ1.w)))

%% (G) The triangular FullJacobian must build the same Jacobian as the plain loop
% FullJacobianReuseVpath=1 restarts the backward pass at the perturbed period rather than redoing all of
% it. That is exact, not an approximation, so the two must agree to machine precision. maxiter=1 is all
% this needs: one Jacobian gets built and one step taken. Run with and without fastOLG, since the
% baseline value fn path it stores has a different shape in each.
fastOLGopts=[1,0];
for ff=1:length(fastOLGopts)
    transpathoptionsG=transpathoptionsGE;
    transpathoptionsG.fastOLG=fastOLGopts(ff);
    transpathoptionsG.maxiter=1;
    transpathoptionsG.GEnewprice=1;
    transpathoptionsG.verbose=0;
    transpathoptionsG.GEnewprice1.Jacobianmethod='FullJacobian';
    transpathoptionsG.GEnewprice1.BroydenRegularisation='minimumnorm';
    transpathoptionsG.GEnewprice1.FullJacobianReuseVpath=0; % the plain loop, one full path solve per column
    tic;
    [PricePathG1,~]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsG, simoptions, vfoptions);
    time_G1=toc;
    transpathoptionsG.GEnewprice1.FullJacobianReuseVpath=1; % restart the backward pass at the perturbed period
    tic;
    [PricePathG2,~]=TransitionPath_Case1_FHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsG, simoptions, vfoptions);
    time_G2=toc;
    fprintf('FullJacobian plain loop vs triangular (fastOLG=%i), one iteration, this should be zero, r: %.3e, w: %.3e \n',fastOLGopts(ff),max(abs(PricePathG2.r-PricePathG1.r)),max(abs(PricePathG2.w-PricePathG1.w)))
    fprintf('Runtime, FullJacobian one iteration (fastOLG=%i): plain loop %2.2f seconds, triangular %2.2f seconds \n',fastOLGopts(ff),time_G1,time_G2)
end

%% (J) The fake-news Jacobian against the brute-force one
% LudwigSSJ builds the sequence-space Jacobian at the final stationary eqm by the fake-news algorithm;
% FullJacobian finite-differences the path itself, one full solve per (period,price). At a path close
% to the stationary eqm the two are derivatives of the same map at nearly the same point, so they must
% nearly agree. Compared through the one quasi-Newton step each takes from the same start: the step is
% -factor*J\f with the same f, so the steps differ only by the difference in J. The start is the bumped
% guess shrunk to 5% of its size, so that the path Jacobian is taken close to the steady state that the
% fake-news one is taken at; what is left is that 5% of nonlinearity plus the finite-difference error of
% both, which is why this is 'small' and not 'zero'. LudwigPath, which keeps none of the cross-time
% structure, is printed alongside as the yardstick for what a genuinely different Jacobian does.
PricePathSmall.r=PricePathGE.r+0.05*(PricePathBumped.r-PricePathGE.r);
PricePathSmall.w=PricePathGE.w+0.05*(PricePathBumped.w-PricePathGE.w);
transpathoptionsJ=transpathoptionsGE;
transpathoptionsJ.GEnewprice=1;
transpathoptionsJ.maxiter=1;
transpathoptionsJ.verbose=0;
transpathoptionsJ.GEnewprice1.BroydenRegularisation='minimumnorm';
transpathoptionsJ.GEnewprice1.Jacobianmethod='FullJacobian';
transpathoptionsJ.GEnewprice1.FullJacobianReuseVpath=1;
PricePathJFull=TransitionPath_Case1_FHorz(PricePathSmall, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsJ, simoptions, vfoptions);
transpathoptionsJ.GEnewprice1.Jacobianmethod='LudwigSSJ';
PricePathJSSJ=TransitionPath_Case1_FHorz(PricePathSmall, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsJ, simoptions, vfoptions);
transpathoptionsJ.GEnewprice1.Jacobianmethod='LudwigPath';
PricePathJPath=TransitionPath_Case1_FHorz(PricePathSmall, ParamPathGE, T, V_finalGE, AgentDist_initialGE, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsJ, simoptions, vfoptions);
stepFull=[PricePathJFull.r-PricePathSmall.r, PricePathJFull.w-PricePathSmall.w];
stepSSJ=[PricePathJSSJ.r-PricePathSmall.r, PricePathJSSJ.w-PricePathSmall.w];
stepPath=[PricePathJPath.r-PricePathSmall.r, PricePathJPath.w-PricePathSmall.w];
fprintf('Size of the FullJacobian quasi-Newton step from the small bump (max abs): %.3e (if this is zero the step was rejected and the comparison below means nothing) \n',max(abs(stepFull)))
fprintf('LudwigSSJ vs FullJacobian, one quasi-Newton step from the small bump, relative difference in the step, should be small: %.3e \n',max(abs(stepSSJ-stepFull))/max(abs(stepFull)))
fprintf('For comparison, LudwigPath vs FullJacobian on the same step, relative difference: %.3e \n',max(abs(stepPath-stepFull))/max(abs(stepFull)))

%% A graph of the solutions, so that a check that fails can be looked at
figure(figure_c);
subplot(2,1,1); plot(1:T,PricePathBumped.r,1:T,PricePathB.r,1:T,PricePathA.r,1:T,PricePathGE.r,'--')
title('r: bumped start, shooting, last quasi-Newton (LudwigSSJ), eqm')
subplot(2,1,2); plot(1:T,PricePathBumped.w,1:T,PricePathB.w,1:T,PricePathA.w,1:T,PricePathGE.w,'--')
title('w: bumped start, shooting, last quasi-Newton (LudwigSSJ), eqm')
legend('start','shooting','quasi-Newton','eqm')
% saved by the driver, as in every Core bank

output=1;

end
