function output=EstimFHorzPTypeMLE_LogLik(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions)
% doPart(1): LogLikelihood_FHorz_PType, the finite mixture over the permanent types, against the single-type LogLikelihood_FHorz called
% by hand type by type (with that type's policy and parameters), and the mixture, posterior and conditional likelihood computed by hand.
%   (a) LogLi_bytype vs single-type calls; LogLi vs log(sum_k w_k exp(LogLi^k)); LogL; the posterior type probabilities
%   (b) masses [1;0]: LogLi equals the single-type likelihood of the first type, posterior one-hot, the zero-mass type all -Inf
%   (c) Data.ptype (observed for most individuals, NaN for every fifth): log(w_k)+LogLi^k where observed, the mixture where not
%   (d) initialcondition='conditional': mixture minus the log of the mixture of the first-age terms
%   (e) how well the posterior picks the true type (the types differ by 50 percent in earnings)
%   (f) jequaloneDist as a structure over the types, and with the type as its last dimension: the same LogLi
%   (g) refusals
% The data: one panel simulated by hand from the model at the truth (EstimFHorzPTypeMLE_SimulatePTypePanel, fixed seed).

Params.sig_c=0.1; Params.sig_y=0.05; Params.sig_a=0.2; % the truth (measurement error)
FnsToEvaluate.cons=@(aprime,a,z,r,w,kappa_j_pt,agej,Jr,pension) (1+r)*a+(agej<Jr)*w*kappa_j_pt*z+(agej>=Jr)*pension-aprime;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt,agej,Jr) (agej<Jr)*w*kappa_j_pt*z;
FnsToEvaluate.assets=@(aprime,a,z) a;
MeasurementError.cons.distribution='normal'; MeasurementError.cons.sigma='sig_c';
MeasurementError.earnings.distribution='normal'; MeasurementError.earnings.sigma='sig_y';
MeasurementError.assets.distribution='normal'; MeasurementError.assets.sigma='sig_a';
N_i=length(Names_i);

[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
[Data,truetype]=EstimFHorzPTypeMLE_SimulatePTypePanel(Policy,n_a,n_z,N_j,Names_i,a_grid,z_grid,pi_z,jequaloneDist,Params,vfoptions.ngridinterp,400,20261016);
N=size(Data.cons,1);
w=Params.(PTypeDistParamNames{1})(:);

% By hand: each type with the single-type command
LogLi_k=zeros(N,N_i);
First_k=zeros(N,N_i);
for ii=1:N_i
    P=Params;
    P.kappa_j_pt=Params.kappa_j_pt(ii,:);
    P.(PTypeDistParamNames{1})=w(ii);
    [~,temp,tempij]=LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,Policy.(Names_i{ii}),FnsToEvaluate,P,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions);
    LogLi_k(:,ii)=gather(temp);
    First_k(:,ii)=gather(tempij(:,1)); % everyone is observed from age 1
end
handmix=log(sum(w'.*exp(LogLi_k-max(LogLi_k,[],2)),2))+max(LogLi_k,[],2);
handpost=w'.*exp(LogLi_k-handmix);

%% (a)
[LogL,LogLi,Post,LogLi_bytype]=LogLikelihood_FHorz_PType(Data,MeasurementError,jequaloneDist,PTypeDistParamNames,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions);
fprintf('PTypeLogLik (a) LogLi_bytype vs single-type calls, this should be zero: %.3e; LogLi vs the mixture by hand, this should be zero: %.3e \n',CoreMaxAbs(LogLi_bytype-LogLi_k),CoreMaxAbs(LogLi-handmix))
fprintf('PTypeLogLik (a) LogL vs sum of LogLi, this should be zero: %.3e; posterior vs by hand, this should be zero: %.3e; posterior rows sum to one, this should be zero: %.3e \n',abs(LogL-sum(handmix)),CoreMaxAbs(Post-handpost),CoreMaxAbs(sum(Post,2)-1))

%% (b) all the mass on the first type
Params_b=Params; Params_b.(PTypeDistParamNames{1})=[1;0];
[~,LogLi_b,Post_b,LogLi_bytype_b]=LogLikelihood_FHorz_PType(Data,MeasurementError,jequaloneDist,PTypeDistParamNames,Policy,FnsToEvaluate,Params_b,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions);
fprintf('PTypeLogLik (b) masses [1;0]: LogLi vs the first type alone, this should be zero: %.3e; posterior of the first type minus one, this should be zero: %.3e; finite entries for the zero-mass type, this should be zero: %i \n',CoreMaxAbs(LogLi_b-LogLi_k(:,1)),CoreMaxAbs(Post_b(:,1)-1),sum(isfinite(LogLi_bytype_b(:,2))))

%% (c) observed types
Data_c=Data;
Data_c.ptype=truetype;
Data_c.ptype(1:5:end)=NaN;
[~,LogLi_c,Post_c]=LogLikelihood_FHorz_PType(Data_c,MeasurementError,jequaloneDist,PTypeDistParamNames,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions);
hand_c=handmix;
known=~isnan(Data_c.ptype);
hand_c(known)=log(w(Data_c.ptype(known)))+LogLi_k(sub2ind([N,N_i],find(known),Data_c.ptype(known)));
onehot=zeros(sum(known),N_i);
onehot(sub2ind([sum(known),N_i],(1:sum(known))',Data_c.ptype(known)))=1;
fprintf('PTypeLogLik (c) observed types: LogLi vs by hand, this should be zero: %.3e; posterior of the observed rows vs one-hot, this should be zero: %.3e; of the unobserved rows vs the mixture posterior, this should be zero: %.3e \n',CoreMaxAbs(LogLi_c-hand_c),CoreMaxAbs(Post_c(known,:)-onehot),CoreMaxAbs(Post_c(~known,:)-handpost(~known,:)))

%% (d) conditional on the first age
estimoptions.initialcondition='conditional';
[~,LogLi_d]=LogLikelihood_FHorz_PType(Data,MeasurementError,jequaloneDist,PTypeDistParamNames,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions,estimoptions);
handfirst=log(sum(w'.*exp(First_k-max(First_k,[],2)),2))+max(First_k,[],2);
fprintf('PTypeLogLik (d) conditional: LogLi vs the mixture minus the mixture of the first-age terms, this should be zero: %.3e \n',CoreMaxAbs(LogLi_d-(handmix-handfirst)))
fprintf('PTypeLogLik (d) (discriminator) the mixture of the first-age terms is as large as %.4f in absolute value (informational) \n',max(abs(handfirst)))

%% (e) the posterior picks the true type
posttrue=Post(sub2ind([N,N_i],(1:N)',truetype));
fprintf('PTypeLogLik (e) one minus the mean posterior probability of the true type, this should be close to zero: %.3e (share classified correctly %.4f, informational) \n',1-mean(posttrue),mean(posttrue>0.5))

%% (f) other forms of jequaloneDist
jeq_f.(Names_i{1})=jequaloneDist; jeq_f.(Names_i{2})=jequaloneDist;
[~,LogLi_f1]=LogLikelihood_FHorz_PType(Data,MeasurementError,jeq_f,PTypeDistParamNames,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions);
jeq_f2=cat(3,w(1)*jequaloneDist,w(2)*jequaloneDist); % the type as the last dimension (its slice masses are the type masses)
simoptions_f=simoptions; simoptions_f.warnjequaloneptypeasdim=0;
[~,LogLi_f2]=LogLikelihood_FHorz_PType(Data,MeasurementError,jeq_f2,PTypeDistParamNames,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_f);
fprintf('PTypeLogLik (f) jequaloneDist as a structure over the types, this should be zero: %.3e; with the type as the last dimension, this should be zero: %.3e \n',CoreMaxAbs(LogLi_f1-LogLi),CoreMaxAbs(LogLi_f2-LogLi))

%% (g) refusals
errcases={'Data.ptype out of range','type masses of the wrong length','negative type masses','bad initialcondition'};
errmsgs={'Data.ptype must be integers from 1 to N_i','does not have the same number of permanent types','must be non-negative','must be ''model'' or ''conditional'''};
for ee=1:length(errcases)
    D=Data; P=Params; esto=struct();
    switch ee
        case 1
            D.ptype=truetype; D.ptype(3)=N_i+1;
        case 2
            P.(PTypeDistParamNames{1})=[0.5;0.3;0.2];
        case 3
            P.(PTypeDistParamNames{1})=[1.2;-0.2];
        case 4
            esto.initialcondition='first';
    end
    errok=0;
    try
        LogLikelihood_FHorz_PType(D,MeasurementError,jequaloneDist,PTypeDistParamNames,Policy,FnsToEvaluate,P,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions,esto);
    catch err
        errok=contains(err.message,errmsgs{ee});
        if ~errok
            fprintf('   (the error was: %s) \n',err.message);
        end
    end
    fprintf('PTypeLogLik (g) %s errors with the intended message, this should be zero: %i \n',errcases{ee},1-errok);
end

output.LogLi=LogLi;
end
