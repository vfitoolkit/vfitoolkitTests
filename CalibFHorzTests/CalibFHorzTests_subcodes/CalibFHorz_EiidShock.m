function output=CalibFHorz_EiidShock(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% CalibFHorzTests part 7: an iid shock e given through vfoptions.EiidShockFn with its dispersion calibrated (the e grids are then
% rebuilt every evaluation), and both shock functions together with rho_z and sigma_e calibrated (either alone must trigger the
% rebuild; before 2026-10-07 EiidShockFn only did when there was no ExogShockFn). Targets: the age-conditional standard deviation
% of earnings while working and the Gini of earnings, generated from a solve with the same e grid at the truth. The grid
% interpolation layer is on. Each run is wrapped so an error in the plumbing is reported as a failed check.
% (ReturnFn and jequaloneDist of the bank's driver are for the model without e and unused here.)

ReturnFn_e=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate_e.assets=@(aprime,a,z,e) a;
FnsToEvaluate_e.earnings=@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;
Jr=Params.Jr;
n_e=3;
Params_3=Params; Params_3.sigma_e=0.1; Params_3.rho_z=0.9; Params_3.sigma_z=0.03;
[e_grid,pi_e]=CalibFHorz_EiidShockFn(Params_3.sigma_e);
vfoptions_e=struct(); vfoptions_e.gridinterplayer=1; vfoptions_e.ngridinterp=15; vfoptions_e.n_e=n_e; vfoptions_e.e_grid=e_grid; vfoptions_e.pi_e=pi_e;
simoptions_e=struct(); simoptions_e.gridinterplayer=1; simoptions_e.ngridinterp=15; simoptions_e.n_e=n_e; simoptions_e.e_grid=e_grid; simoptions_e.pi_e=pi_e;
simoptions_e_full=simoptions_e; simoptions_e_full.whichstats=ones(1,7);
jequaloneDist_e=zeros(n_a,n_z,n_e,'gpuArray'); jequaloneDist_e(1,ceil(n_z/2),ceil(n_e/2))=1;
caliboptions=struct();
caliboptions.verbose=0;
caliboptions.fminalgo=8;

ran=0;
try
    [~,Policy_3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_3,DiscountFactorParamNames,[],vfoptions_e);
    StationaryDist_3=StationaryDist_FHorz_Case1(jequaloneDist_e,AgeWeightParamNames,Policy_3,n_d,n_a,n_z,N_j,pi_z,Params_3,simoptions_e);
    AllFull_3=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_3,Policy_3,FnsToEvaluate_e,Params_3,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_e_full);
    LCPFull_3=LifeCycleProfiles_FHorz_Case1(StationaryDist_3,Policy_3,FnsToEvaluate_e,Params_3,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_e_full);
    ran=1;
catch err
    fprintf('EiidShock: the model solve or the stats commands errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('EiidShock: the model with e solved and its stats ran, did not run, this should be zero: %.3e \n',1-ran)
if ran==1
    TM3=struct();
    t=gather(LCPFull_3.earnings.StdDeviation); t(Jr:end)=NaN; TM3.AgeConditionalStats.earnings.StdDeviation=t;
    TM3.AllStats.earnings.Gini=gather(AllFull_3.earnings.Gini);
    vfoptions_c=struct(); vfoptions_c.gridinterplayer=1; vfoptions_c.ngridinterp=15; vfoptions_c.n_e=n_e; vfoptions_c.EiidShockFn=@(sigma_e) CalibFHorz_EiidShockFn(sigma_e);
    simoptions_c=struct(); simoptions_c.gridinterplayer=1; simoptions_c.ngridinterp=15; simoptions_c.n_e=n_e; simoptions_c.EiidShockFn=vfoptions_c.EiidShockFn;
    for ss=1:2
        Params_start=Params_3; lab='sigma_e started at the truth';
        if ss==2
            Params_start.sigma_e=0.15; lab='sigma_e started at 0.15 (truth 0.1)';
        end
        ran=0; objv=NaN; perr=NaN; tt=NaN;
        try
            t0=tic;
            [CalibParams,calibsummary]=CalibrateLifeCycleModel({'sigma_e'},TM3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_start,DiscountFactorParamNames,jequaloneDist_e,AgeWeightParamNames,[],FnsToEvaluate_e,caliboptions,vfoptions_c,simoptions_c);
            tt=toc(t0);
            ran=1; objv=calibsummary.objvalue; perr=abs(CalibParams.sigma_e-0.1);
        catch err
            fprintf('EiidShock, %s: the calibration errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
        end
        if ss==1
            fprintf('EiidShock, EiidShockFn, %s (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; sigma_e moved, this should be zero: %.3e \n',lab,tt,1-ran,objv,perr)
        else
            fprintf('EiidShock, EiidShockFn, %s (%.1f s), objective %.3e: did not run, this should be zero: %.3e; sigma_e error, this should be close to zero: %.3e \n',lab,tt,objv,1-ran,perr)
        end
    end
    vfoptions_c.ExogShockFn=@(rho_z,sigma_z) CalibFHorz_ExogShockFn(rho_z,sigma_z);
    simoptions_c.ExogShockFn=vfoptions_c.ExogShockFn;
    Params_start=Params_3; Params_start.rho_z=0.8; Params_start.sigma_e=0.15;
    caliboptions_both=caliboptions; caliboptions_both.constrain0to1={'rho_z'}; % keep the autocorrelation in (0,1), as a user would: lsqnonlin's finite-difference steps took it past one
    ran=0; objv=NaN; perr=NaN; tt=NaN;
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateLifeCycleModel({'rho_z','sigma_e'},TM3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_start,DiscountFactorParamNames,jequaloneDist_e,AgeWeightParamNames,[],FnsToEvaluate_e,caliboptions_both,vfoptions_c,simoptions_c);
        tt=toc(t0);
        ran=1; objv=calibsummary.objvalue; perr=max(abs([CalibParams.rho_z;CalibParams.sigma_e]-[0.9;0.1]));
    catch err
        fprintf('EiidShock, ExogShockFn and EiidShockFn together: the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('EiidShock, ExogShockFn and EiidShockFn together, rho_z and sigma_e from (0.8, 0.15) (truth 0.9, 0.1) (%.1f s), objective %.3e: did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',tt,objv,1-ran,perr)
end

output=struct();
output.ran=ran;

end
