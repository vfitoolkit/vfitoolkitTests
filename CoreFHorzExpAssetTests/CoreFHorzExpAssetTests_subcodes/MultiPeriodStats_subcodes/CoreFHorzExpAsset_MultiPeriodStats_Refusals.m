function output=CoreFHorzExpAsset_MultiPeriodStats_Refusals(n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline)
% The refusal set of the multi-period statistics: EvalFnOnAgentDist_AutoCorrTransProbs_FHorz, and the M>1 windows of
% EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz, EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz and
% LifeCycleProfiles_FHorz_Case1, must each refuse (with their 'not yet implemented' error) every asset whose next-period value
% is not read from the Policy and that they do not yet handle: experienceassetu, experienceassetz, experienceassete,
% experienceassetze, experienceassetsemiz, riskyasset, residualasset, inheritanceasset. Checked as a set (every command with
% every flag), because a refusal that is missing for one of them is the case that silently gives wrong numbers (before these
% refusals, AutoCorrTransProbs with gridinterplayer read the L2 index row of the Policy as a2').
% The experienceasset itself is handled by all four (CoreFHorzExpAsset_MultiPeriodStats_FakeExpAsset and _BruteForce check
% the numbers), so here it is the positive control: the same calls with experienceasset=1 must run.
% The model is the fake experience asset (aprimeFn=@(d2,a2) d2) with z and e; each flag is put on its simoptions in place of
% experienceasset (all of these assets are the last endogenous state, so the Policy has the same shape), except
% inheritanceasset, which is put on the standard two-endogenous-state version of the same model.

fprintf('\n========== CoreFHorzExpAsset_MultiPeriodStats_Refusals ==========\n')

aprimeFn=@(d2,a2) d2;
d_grid=a_grid(n_a(1)+1:end);
n_d=n_a(2);

vfoptions0=struct();
vfoptions0.n_e=vfoptionsbaseline.n_e;
vfoptions0.e_grid=vfoptionsbaseline.e_grid;
vfoptions0.pi_e=vfoptionsbaseline.pi_e;
simoptions0=vfoptions0;
simoptions0.pi_z=pi_z;
ReturnFn0=@(a1prime,a2prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_TwoEndo_nod1_noz_noe_nosemiz(a1prime,a2prime,a1,a2,r,w*z*e,kappa_j,sigma,agej,Jr,pension);
ReturnFn1=@(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz(d2,a1prime,a1,a2,r,w*z*e,kappa_j,sigma,agej,Jr,pension);
F0.assets=@(a1prime,a2prime,a1,a2,z,e) a1;
F0.human=@(a1prime,a2prime,a1,a2,z,e) a2;
F1.assets=@(d2,a1prime,a1,a2,z,e) a1;
F1.human=@(d2,a1prime,a1,a2,z,e) a2;
jequaloneDist=zeros([n_a,n_z,vfoptions0.n_e],'gpuArray');
jequaloneDist(1,ceil(n_a(2)/2),ceil(n_z/2),ceil(vfoptions0.n_e/2))=1; % a2 mid-grid (at a2=0 every choice is -Inf)

vfoptions1=vfoptions0;
vfoptions1.experienceasset=1;
vfoptions1.aprimeFn=aprimeFn;
simoptions1=simoptions0;
simoptions1.experienceasset=1;
simoptions1.aprimeFn=aprimeFn;
simoptions1.d_grid=d_grid;
simoptions1.a_grid=a_grid;

[~,Policy0]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn0,Params,DiscountFactorParamNames,[],vfoptions0);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy0,0,n_a,n_z,N_j,pi_z,Params,simoptions0);
[~,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn1,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptions1);

flags={'experienceassetu','experienceassetz','experienceassete','experienceassetze','experienceassetsemiz','riskyasset','residualasset','inheritanceasset'};
cmds={'AutoCorrTransProbs_FHorz','CrossSectionCovarCorr_FHorz (M=2)','AgeConditionalStats_CrossSectionCovarCorr_FHorz (M=2)','LifeCycleProfiles_FHorz_Case1 (M=2)'};
nchecks=0;
for kk=0:length(flags)
    for cmd=1:length(cmds)
        if kk==0 % the positive control: the experience asset itself runs
            so=simoptions1; Pol=Policy1; SD=StationaryDist1; F=F1; nd=n_d; dg=d_grid; flagstr='experienceasset';
        elseif strcmp(flags{kk},'inheritanceasset')
            so=simoptions0; Pol=Policy0; SD=StationaryDist0; F=F0; nd=0; dg=[]; flagstr=flags{kk};
            so.(flags{kk})=1;
        else
            so=rmfield(simoptions1,'experienceasset'); Pol=Policy1; SD=StationaryDist1; F=F1; nd=n_d; dg=d_grid; flagstr=flags{kk};
            so.(flags{kk})=1;
        end
        so.ntimeperiods=[1,2];
        errmsg='';
        try
            if cmd==1
                EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(SD,Pol,F,Params,[],nd,n_a,n_z,N_j,dg,a_grid,z_grid,pi_z,so);
            elseif cmd==2
                EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(SD,Pol,F,Params,[],nd,n_a,n_z,N_j,dg,a_grid,z_grid,so);
            elseif cmd==3
                EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(SD,Pol,F,Params,[],nd,n_a,n_z,N_j,dg,a_grid,z_grid,so);
            else
                LifeCycleProfiles_FHorz_Case1(SD,Pol,F,Params,[],nd,n_a,n_z,N_j,dg,a_grid,z_grid,so);
            end
        catch ME
            errmsg=ME.message;
        end
        nchecks=nchecks+1;
        if kk==0
            fprintf('MultiPeriodStats refusal set: %s with experienceasset runs (positive control), errored, this should be zero: %i \n',cmds{cmd},~isempty(errmsg))
        else
            refused=contains(errmsg,'not yet implemented') && contains(errmsg,flagstr);
            fprintf('MultiPeriodStats refusal set: %s with %s refused, not refused, this should be zero: %i \n',cmds{cmd},flagstr,~refused)
        end
        if ~isempty(errmsg) && (kk==0 || ~refused)
            fprintf('    the error was: %s \n',errmsg)
        end
    end
end
fprintf('MultiPeriodStats refusal set: %i (command, flag) pairs checked, of %i \n',nchecks,length(cmds)*(1+length(flags)))

output=struct();

end
