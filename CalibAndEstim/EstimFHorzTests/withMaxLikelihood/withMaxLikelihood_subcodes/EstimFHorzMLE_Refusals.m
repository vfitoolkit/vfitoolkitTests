function output=EstimFHorzMLE_Refusals()
% doPart(3): inputs LogLikelihood_FHorz must refuse, each with its intended message.
% A tiny hand-built model (one endogenous state, no shocks, N_j=4); every case first runs the valid baseline, then changes one thing.
% Each case prints 1-errok, so a case that does NOT error (or errors with a different message) fails.

N_j=4;
n_d=0; d_grid=[];
n_a=6; a_grid=linspace(0,2,n_a)';
n_z=0; z_grid=[]; pi_z=[];
Params.r=0.05; Params.w=1; Params.kappa_j=[1,1.2,1.4,0.8];
Params.sig_c=0.3; Params.sig_neg=-0.3; Params.pmis=0.1;
Policy=gpuArray(repmat(reshape(mod((1:n_a)',n_a)+1,[1,n_a]),1,1,N_j)); % aprime index: a+1, wrapping to 1
jequaloneDist=zeros(n_a,1); jequaloneDist(1)=1;
FnsToEvaluate.cons=@(aprime,a,r,w,kappa_j) (1+r)*a+w*kappa_j-aprime;
FnsToEvaluate.aprimeobs=@(aprime,a) aprime;
simoptions=struct();

Data0.cons=[1,1.1,1.2; 0.9,NaN,1.3];
Data0.age=[1,2,3; 2,3,4];
ME0.cons.distribution='normal'; ME0.cons.sigma='sig_c';

% the valid baseline runs
LogL=LogLikelihood_FHorz(Data0,ME0,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions);
fprintf('Refusals: the valid baseline runs, LogL=%.4f (informational) \n',gather(LogL))

errcases={'no ages','age out of range','two observations at one age','observation where age is NaN','agej0 runs past N_j with data', ...
    'no measurement error for a variable','measurement error for a variable not in Data','sigma not in Params','negative sigma', ...
    'log with non-positive data','log with exact','unknown distribution','misclassified datum not a possible value', ...
    'misclassification matrix without misclassmatrix','experienceasset not implemented','no FnsToEvaluate for a variable', ...
    'bad initialcondition','jequaloneDist not of mass one'};
errmsgs={'must contain the ages','must be integers from 1 to N_j','two observations at the same age','has observations in periods where the age is NaN', ...
    'has observations in periods where the age is NaN','needs a measurement error specification','has no corresponding observed variable in Data', ...
    'which is not in Parameters','must be positive','non-positive observations','only possible with the normal and t distributions', ...
    'distribution must be one of','not among the possible values','misclassmatrix is required','is not yet implemented', ...
    'has no corresponding FnsToEvaluate','must be ''model'' or ''conditional''','must be of mass one'};
for ee=1:length(errcases)
    Data=Data0; ME=ME0; simo=simoptions; esto=struct(); jeq=jequaloneDist; Fns=FnsToEvaluate;
    switch ee
        case 1
            Data=rmfield(Data,'age');
        case 2
            Data.age(2,3)=N_j+1;
        case 3
            Data.age(1,3)=2;
        case 4
            Data.age(1,2)=NaN;
        case 5
            Data=rmfield(Data,'age'); Data.agej0=[1;3]; % individual 2: ages 3,4,5, and 5>N_j has an observation
        case 6
            Data.aprimeobs=[0.4,0.6,0.8; NaN,NaN,NaN];
        case 7
            ME.aprimeobs.distribution='exact';
        case 8
            ME.cons.sigma='sig_notthere';
        case 9
            ME.cons.sigma='sig_neg';
        case 10
            ME.cons.type='log'; Data.cons(1,1)=-0.1;
        case 11
            ME.cons.type='log'; ME.cons.distribution='exact';
        case 12
            ME.cons.distribution='laplace';
        case 13
            Data.aprimeobs=[0.4,0.8,1.2; NaN,NaN,0.55]; ME.aprimeobs.distribution='misclassified'; ME.aprimeobs.misclassprob='pmis';
        case 14
            Data.aprimeobs=[0.4,0.8,1.2; NaN,NaN,NaN]; ME.aprimeobs.distribution='misclassified'; ME.aprimeobs.misclassprob='pmis'; ME.aprimeobs.misclassification='matrix';
        case 15
            simo.experienceasset=1;
        case 16
            Fns=rmfield(Fns,'cons');
        case 17
            esto.initialcondition='first';
        case 18
            jeq=2*jequaloneDist;
    end
    errok=0;
    try
        LogLikelihood_FHorz(Data,ME,jeq,Policy,Fns,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simo,esto);
    catch err
        errok=contains(err.message,errmsgs{ee});
        if ~errok
            fprintf('   (the error was: %s) \n',err.message);
        end
    end
    fprintf('Refusals (%i) %s errors with the intended message, this should be zero: %i \n',ee,errcases{ee},1-errok);
end

output.nrefusals=length(errcases);
end
