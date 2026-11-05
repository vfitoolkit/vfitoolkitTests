function output=UncertainTPathFHorz_TreeParams(T,N_j)
% Stage F0: UncertainTransPathTree with parameters that are vectors (e.g. age-dependent, N_j-by-T) and the per-S form
% ParamPath.(name).S1, .S2, ... Written 2026-11-05, before the first run. No solves.
%   (a) the per-S form of a scalar is the N_S-by-T form, exactly
%   (b) a per-S vector (N_j-by-T fields): every node's column is the field of its S at its date, by hand
%   (c) a vector that does not depend on S (N_j-by-T, and T-by-N_j): every node's column is that date's, by hand
%   (d) on a single history, the per-S vector is the deterministic path of the history's S, column by column
%   (e) refusals, each by its message (vectors without the internal flag, a missing S field, fields of different sizes)
% transpathoptions.uncertaintreevectorparams=1 is internal (the FHorz commands set it); it is set here by hand.

tpV.uncertaintreevectorparams=1;
pi_S=[0.75,0.25;0,1]; T1=2; T2=5;
tvec=1:1:T;
jvec=(1:1:N_j)';

%% (a) The per-S form of a scalar
xS=[0.5+0.01*sin(tvec/3); 0.4+0.02*cos(tvec/5)];
ParamPathA.x=xS;
ParamPathB.x.S1=xS(1,:); ParamPathB.x.S2=xS(2,:);
[Tree,~,NodesA]=UncertainTransPathTree(T,T1,T2,pi_S,1,[],ParamPathA);
[~,~,NodesB]=UncertainTransPathTree(T,T1,T2,pi_S,1,[],ParamPathB);
fprintf('Uncertain FHorz tree (a): the per-S form of a scalar is the N_S-by-T form, this should be zero: %.3e \n',max(abs(NodesA.x-NodesB.x)))
fprintf('Uncertain FHorz tree (a): and is 1-by-nNodes without the internal flag, should be one: %i \n',isequal(size(NodesB.x),[1,Tree.nNodes]))

%% (b) A per-S vector, N_j-by-T fields (e.g. a pension that a repeal cuts, by age)
pen1=(jvec>=10).*(0.3+0.001*tvec); % N_j-by-T
pen2=(jvec>=10).*(0.2+0.002*tvec);
ParamPathC.pension.S1=pen1; ParamPathC.pension.S2=pen2;
[~,~,NodesC]=UncertainTransPathTree(T,T1,T2,pi_S,1,[],ParamPathC,tpV);
handC=zeros(N_j,Tree.nNodes);
for kk=1:Tree.nNodes
    if Tree.nodeS(kk)==1
        handC(:,kk)=pen1(:,Tree.nodedate(kk));
    else
        handC(:,kk)=pen2(:,Tree.nodedate(kk));
    end
end
fprintf('Uncertain FHorz tree (b): a per-S N_j-by-T parameter is N_j-by-nNodes, should be one: %i \n',isequal(size(NodesC.pension),[N_j,Tree.nNodes]))
fprintf('Uncertain FHorz tree (b): every node''s column is the field of its S at its date, this should be zero: %.3e \n',max(abs(NodesC.pension-handC),[],'all'))
fprintf('Uncertain FHorz tree (b): for scale, the two fields differ by, should NOT be zero: %.3e \n',max(abs(pen1-pen2),[],'all'))

%% (c) A vector that does not depend on S
agew=repmat(jvec/sum(jvec),1,T).*(1+0.001*tvec); % N_j-by-T
ParamPathD.agew=agew;
[~,~,NodesD]=UncertainTransPathTree(T,T1,T2,pi_S,1,[],ParamPathD,tpV);
ParamPathD2.agew=agew';
[~,~,NodesD2]=UncertainTransPathTree(T,T1,T2,pi_S,1,[],ParamPathD2,tpV);
fprintf('Uncertain FHorz tree (c): an N_j-by-T parameter, every node''s column is that date''s, this should be zero: %.3e \n',max(abs(NodesD.agew-agew(:,Tree.nodedate)),[],'all'))
fprintf('Uncertain FHorz tree (c): and T-by-N_j is read the same, this should be zero: %.3e \n',max(abs(NodesD2.agew-NodesD.agew),[],'all'))

%% (d) A single history: the per-S vector is the path of the history's S
pi_S1=[0,1;0,1];
[Tree1,~,Nodes1]=UncertainTransPathTree(T,1,1,pi_S1,1,[],ParamPathC,tpV);
Spath=Tree1.historyS(1,:);
hand1=zeros(N_j,T);
hand1(:,Spath==1)=pen1(:,Spath==1);
hand1(:,Spath==2)=pen2(:,Spath==2);
fprintf('Uncertain FHorz tree (d): a single history (S=1 then 2), the per-S vector is the history''s path, this should be zero: %.3e \n',max(abs(Nodes1.pension-hand1),[],'all'))

%% (e) Refusals, each by its message
err=zeros(1,4);
try
    UncertainTransPathTree(T,T1,T2,pi_S,1,[],ParamPathC);
catch ME
    err(1)=contains(ME.message,'only the FHorz uncertain transition path commands take parameters that are vectors');
end
try
    UncertainTransPathTree(T,T1,T2,pi_S,1,[],ParamPathD);
catch ME
    err(2)=contains(ME.message,'or given per S as ParamPath.agew.S1');
end
ParamPathE.pension.S1=pen1;
try
    UncertainTransPathTree(T,T1,T2,pi_S,1,[],ParamPathE,tpV);
catch ME
    err(3)=contains(ME.message,'is given per S, but has no field S2');
end
ParamPathE.pension.S2=pen2(1:N_j-1,:);
try
    UncertainTransPathTree(T,T1,T2,pi_S,1,[],ParamPathE,tpV);
catch ME
    err(4)=contains(ME.message,'are not all the same size');
end
fprintf('Uncertain FHorz tree (e): a per-S vector without the internal flag (the InfHorz commands) is refused by name, should be one: %i \n',err(1))
fprintf('Uncertain FHorz tree (e): an N_j-by-T parameter without the flag is refused by name, should be one: %i \n',err(2))
fprintf('Uncertain FHorz tree (e): a per-S parameter missing an S that a node is in is refused by name, should be one: %i \n',err(3))
fprintf('Uncertain FHorz tree (e): per-S fields of different sizes are refused by name, should be one: %i \n',err(4))

output=struct();

end
