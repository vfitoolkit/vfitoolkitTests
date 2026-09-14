function output=CoreInfHorzInheritAsset_CrossTests_guards(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline)
% Probes of the model shapes that InfHorz inheritanceasset does NOT support.
%
% ValueFnIter_InfHorz_InheritAsset is a cascade of eight 'Have not yet implemented' errors around
% a single live branch. Each probe below drives one of those branches and records whether it
% errored. Everything is wrapped in try/catch for two reasons: a guard that stops firing must
% show up as a failed check rather than as a silently-accepted bad input, and running these bare
% would abort the bank.
%
% If support is ever added for one of these shapes, the corresponding line flips to 0 and the
% diary says so, which is the signal to come back and write a real test for it.

n_d1=n_d(1);
n_d2=n_d(2);
N_d1=prod(n_d1);
d1_grid=d_grid(1:N_d1);
d2_grid=d_grid(N_d1+1:end);

ReturnFn=@(d1,d2,a,z,r,w,sigma,eta,varphi) ReturnFn_d1d2_a_z(d1,d2,a,z,r,w,sigma,eta,varphi);

fprintf('\n===== inheritanceasset shape guards ===== \n')

%% 0. Positive control: the supported shape solves and gives a finite V
% Without this, every 'errors as expected' below would also be consistent with the whole family
% being broken.
ok=0; msg='';
try
    vf=vfoptionsbaseline;
    Vt=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vf);
    Vt=gather(Vt);
    ok=any(isfinite(Vt(:))) && ~any(isnan(Vt(:)));
catch ME
    msg=ME.message;
end
fprintf('0. supported shape (d1, no a1, z, no e) solves with no NaN, this should be one: %d \n',ok)
if ok~=1
    fprintf('   unexpected error: %s \n',msg)
end

%% 1. No d1 must error
% n_d scalar means n_d1=0, which the dispatcher rejects.
ok=0; msg='';
try
    vf=vfoptionsbaseline;
    ReturnFn_nod1=@(d2,a,z,r,w,sigma,eta,varphi) ReturnFn_d1d2_a_z(1,d2,a,z,r,w,sigma,eta,varphi);
    Vt=ValueFnIter_InfHorz(n_d2,n_a,n_z,d2_grid,a_grid,z_grid,pi_z,ReturnFn_nod1,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=1; msg=ME.message;
end
fprintf('1. no d1 errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 2. With a1 must error
% Two endogenous states means n_a1>0, which the dispatcher rejects.
ok=0; msg='';
try
    vf=vfoptionsbaseline;
    ReturnFn_a1=@(d1,d2,a1prime,a1,a2,z,r,w,sigma,eta,varphi) ReturnFn_d1d2_a_z(d1,d2,a2,z,r,w,sigma,eta,varphi);
    Vt=ValueFnIter_InfHorz(n_d,[11,n_a],n_z,d_grid,[linspace(0,1,11)';a_grid],z_grid,pi_z,ReturnFn_a1,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=1; msg=ME.message;
end
fprintf('2. with a1 errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 3. With e must error
ok=0; msg='';
try
    vf=vfoptionsbaseline;
    vf.n_e=3;
    vf.e_grid=linspace(0.9,1.1,3)';
    vf.pi_e=[1/3;1/3;1/3];
    ReturnFn_e=@(d1,d2,a,z,e,r,w,sigma,eta,varphi) ReturnFn_d1d2_a_z(d1,d2,a,z*e,r,w,sigma,eta,varphi);
    Vt=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=1; msg=ME.message;
end
fprintf('3. with e errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 4. No z must error
% The inheritance asset is defined off the shock transition, so there is nothing for it to do
% without z; the dispatcher says so explicitly.
ok=0; msg='';
try
    vf=vfoptionsbaseline;
    ReturnFn_noz=@(d1,d2,a,r,w,sigma,eta,varphi) ReturnFn_d1d2_a_z(d1,d2,a,1,r,w,sigma,eta,varphi);
    Vt=ValueFnIter_InfHorz(n_d,n_a,0,d_grid,a_grid,[],[],ReturnFn_noz,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=1; msg=ME.message;
end
fprintf('4. no z errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 5. Missing aprimeFn must error
ok=0; msg='';
try
    vf=struct();
    vf.inheritanceasset=1; % but no vf.aprimeFn
    Vt=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=1; msg=ME.message;
end
fprintf('5. missing aprimeFn errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 6. Missing simoptions.aprimeFn on the distribution side must error
ok=0; msg='';
try
    vf=vfoptionsbaseline;
    [~,Policyt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vf);
    so=struct(); so.inheritanceasset=1; % but no so.aprimeFn
    StationaryDistt=StationaryDist_InfHorz(Policyt,n_d,n_a,n_z,pi_z,so,Params,[]); %#ok<NASGU>
catch ME
    ok=1; msg=ME.message;
end
fprintf('6. missing simoptions.aprimeFn errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%%
output=struct();

end
