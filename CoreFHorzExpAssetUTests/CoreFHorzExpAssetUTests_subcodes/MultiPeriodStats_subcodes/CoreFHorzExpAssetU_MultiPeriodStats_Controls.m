function output=CoreFHorzExpAssetU_MultiPeriodStats_Controls(n_d2,d2_grid,n_a,a_grid,n_z,N_j,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline)
% The option checks of the multi-period statistics for experienceassetu, in EvalFnOnAgentDist_AutoCorrTransProbs_FHorz and the
% M>1 windows of EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz, EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz
% and LifeCycleProfiles_FHorz_Case1, checked as a set (every command):
%   positive control: experienceassetu (with n_u, u_grid, pi_u) runs;
%   experienceasset and experienceassetu both set: refused ('only one of');
%   experienceassetu without n_u (or u_grid, or pi_u): refused (says n_u, u_grid and pi_u are needed).
% (That the other unhandled assets are refused is checked in CoreFHorzExpAssetTests doPart(20).)
% The model is the baseline experienceassetu one (no d1, a1, z and e), on the bank's grids.

fprintf('\n========== CoreFHorzExpAssetU_MultiPeriodStats_Controls ==========\n')

n_d=n_d2;
d_grid=d2_grid;
vfoptions=struct();
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
vfoptions.experienceassetu=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions.n_u=vfoptionsbaseline.n_u;
vfoptions.u_grid=vfoptionsbaseline.u_grid;
vfoptions.pi_u=vfoptionsbaseline.pi_u;
simoptions=vfoptions;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;
simoptions.pi_z=pi_z;
ReturnFn=@(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_e_nosemiz(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
F.assets=@(d2,a1prime,a1,a2,z,e) a1;
F.human=@(d2,a1prime,a1,a2,z,e) a2;
jequaloneDist=zeros([n_a,n_z,vfoptions.n_e],'gpuArray');
jequaloneDist(1,ceil(n_a(2)/2),ceil(n_z/2),ceil(vfoptions.n_e/2))=1;

[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

cases={'experienceassetu runs (positive control)','experienceasset and experienceassetu both set','experienceassetu without n_u','experienceassetu without u_grid','experienceassetu without pi_u'};
expected={'','only one of','n_u, simoptions.u_grid and simoptions.pi_u','n_u, simoptions.u_grid and simoptions.pi_u','n_u, simoptions.u_grid and simoptions.pi_u'};
cmds={'AutoCorrTransProbs_FHorz','CrossSectionCovarCorr_FHorz (M=2)','AgeConditionalStats_CrossSectionCovarCorr_FHorz (M=2)','LifeCycleProfiles_FHorz_Case1 (M=2)'};
nchecks=0;
for kk=1:length(cases)
    for cmd=1:length(cmds)
        so=simoptions;
        so.ntimeperiods=[1,2];
        if kk==2
            so.experienceasset=1;
        elseif kk==3
            so=rmfield(so,'n_u');
        elseif kk==4
            so=rmfield(so,'u_grid');
        elseif kk==5
            so=rmfield(so,'pi_u');
        end
        errmsg='';
        try
            if cmd==1
                EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,F,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,so);
            elseif cmd==2
                EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,F,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,so);
            elseif cmd==3
                EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,F,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,so);
            else
                LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,F,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,so);
            end
        catch ME
            errmsg=ME.message;
        end
        nchecks=nchecks+1;
        if kk==1
            bad=~isempty(errmsg);
            fprintf('MultiPeriodStats U controls: %s, %s, errored, this should be zero: %i \n',cmds{cmd},cases{kk},bad)
        else
            bad=~contains(errmsg,expected{kk});
            fprintf('MultiPeriodStats U controls: %s, %s, not refused with the expected error, this should be zero: %i \n',cmds{cmd},cases{kk},bad)
        end
        if bad && ~isempty(errmsg)
            fprintf('    the error was: %s \n',errmsg)
        end
    end
end
fprintf('MultiPeriodStats U controls: %i (command, case) pairs checked, of %i \n',nchecks,length(cmds)*length(cases))

output=struct();

end
