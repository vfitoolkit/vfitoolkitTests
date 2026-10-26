function [Data,truetype]=EstimFHorzPTypeMLE_SimulatePTypePanel(Policy,n_a,n_z,N_j,Names_i,a_grid,z_grid,pi_z,jequaloneDist,Params,ngridinterp,N,seed)
% A hand-coded simulation of one panel from the EstimFHorzPTypeTests model (two permanent types differing in kappa_j_pt), for the
% maximum likelihood tests. round(N*ptypeweights(k)) individuals of each type k, each simulated by EstimFHorzMLE_SimulatePanel (the
% single-type simulator of the EstimFHorzTests/withMaxLikelihood bank, a plain loop with a fixed seed) with that type's policy and
% earnings profile. Params.sig_c can be a structure over the types (a per-type measurement error s.d.).
% (Test-bank only: the toolkit's estimation never simulates.)

N_i=length(Names_i);
Data=struct();
truetype=[];
for ii=1:N_i
    P=Params;
    P.kappa_j=Params.kappa_j_pt(ii,:);
    if isstruct(Params.sig_c)
        P.sig_c=Params.sig_c.(Names_i{ii});
    end
    N_k=round(N*Params.ptypeweights(ii));
    Data_k=EstimFHorzMLE_SimulatePanel(Policy.(Names_i{ii}),n_a,n_z,N_j,a_grid,z_grid,pi_z,jequaloneDist,P,ngridinterp,N_k,seed+ii);
    if ii==1
        Data=Data_k;
    else
        names=fieldnames(Data_k);
        for ff=1:length(names)
            Data.(names{ff})=[Data.(names{ff}); Data_k.(names{ff})];
        end
    end
    truetype=[truetype; ii*ones(N_k,1)];
end

end
