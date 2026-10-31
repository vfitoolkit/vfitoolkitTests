function output=UncertainTPath_TreeStructure(T)
% The tree of histories, checked against counts and probabilities worked out by hand. No solves.

%% (a) n_S=1: one history, numbered 1,...,T, one segment
Tree=UncertainTransPathTree(T,3,5,1,1);
fprintf('Tree (a): n_S=1 has one history, should be one: %i \n',Tree.nHist==1)
fprintf('Tree (a): and T nodes numbered 1,...,T, should be one: %i \n',Tree.nNodes==T && isequal(Tree.histories,1:T) && isequal(Tree.nodedate(:)',1:T))
fprintf('Tree (a): and a single segment, should be one: %i \n',Tree.nSeg==1 && isequal(Tree.segnodes{1},1:T))

%% (b) absorbing S (pi_S=[P,1-P;0,1]): one history that never switches, plus one per switch date
P=0.75; T1=3; T2=6;
pi_S=[P,1-P;0,1];
Tree=UncertainTransPathTree(T,T1,T2,pi_S,1);
nNodesHand=T+sum(T-(T1:T2)+1); % the never-switch path, plus each branch from its switch date s to T
fprintf('Tree (b): absorbing S has T2-T1+2 histories, should be one: %i \n',Tree.nHist==T2-T1+2)
fprintf('Tree (b): and T+sum_s(T-s+1) nodes, should be one: %i \n',Tree.nNodes==nNodesHand)
hh0=find(all(Tree.historyS==1,2));
fprintf('Tree (b): the never-switch history has probability P^(T2-T1+1), this should be zero: %.3e \n',abs(Tree.historyprob(hh0)-P^(T2-T1+1)))
probhand=zeros(T2-T1+1,1); probtree=zeros(T2-T1+1,1);
for s=T1:T2
    hh=find(Tree.historyS(:,s)==2 & (s==1 | Tree.historyS(:,max(s-1,1))==1));
    probtree(s-T1+1)=Tree.historyprob(hh);
    probhand(s-T1+1)=P^(s-T1)*(1-P);
end
fprintf('Tree (b): the history switching at s has probability P^(s-T1)*(1-P), this should be zero: %.3e \n',max(abs(probtree-probhand)))
% Segments: periods 1 to T1-1 (ending at a branching node), one for each never-switch node in T1 to T2-1 (each a
% branching node on its own), the never-switch tail from T2 to T, and one per branch
fprintf('Tree (b): 1+(T2-T1)+1+(T2-T1+1) segments, should be one: %i \n',Tree.nSeg==1+(T2-T1)+1+(T2-T1+1))

%% (c) full support, n_S=2, a window of 3: 8 histories
T1=4; T2=6;
pi_S=[0.5,0.5;0.25,0.75];
Tree=UncertainTransPathTree(T,T1,T2,pi_S,1);
fprintf('Tree (c): full support over a window of 3 has 2^3 histories, should be one: %i \n',Tree.nHist==8)
fprintf('Tree (c): and (T1-1)+2+4+8+8*(T-T2) nodes, should be one: %i \n',Tree.nNodes==(T1-1)+2+4+8+8*(T-T2))
fprintf('Tree (c): the history probabilities sum to one, this should be zero: %.3e \n',abs(sum(Tree.historyprob)-1))
haspar=(Tree.nodeparent>0);
fprintf('Tree (c): every node probability is its parent''s times pi_S, this should be zero: %.3e \n',max(abs(Tree.nodeprob(haspar)-Tree.nodeprob(Tree.nodeparent(haspar)).*Tree.nodecondprob(haspar))))
pihand=pi_S(sub2ind([2,2],Tree.nodeS(Tree.nodeparent(haspar)),Tree.nodeS(haspar)));
inwindow=(Tree.nodedate(haspar)>=T1 & Tree.nodedate(haspar)<=T2);
fprintf('Tree (c): the conditional probability is pi_S(S_parent,S) in the window and one outside, this should be zero: %.3e \n',max(abs(Tree.nodecondprob(haspar)-(inwindow.*pihand+(1-inwindow)))))
% the probabilities of the nodes in each period sum to one
sumbydate=accumarray(Tree.nodedate,Tree.nodeprob);
fprintf('Tree (c): the node probabilities in each period sum to one, this should be zero: %.3e \n',max(abs(sumbydate-1)))
fprintf('Tree (c): nodes are numbered date major, should be one: %i \n',all(diff(Tree.nodedate)>=0))
fprintf('Tree (c): S is S_initial before T1, should be one: %i \n',all(Tree.historyS(:,1:T1-1)==1,'all'))
fprintf('Tree (c): S is frozen after T2, should be one: %i \n',all(Tree.historyS(:,T2+1:T)==Tree.historyS(:,T2),'all'))
fprintf('Tree (c): the 8 histories are the 8 distinct S sequences on the window, should be one: %i \n',size(unique(Tree.historyS(:,T1:T2),'rows'),1)==8)
% every node is in exactly one segment, a segment is a chain, and only its last node can branch
covered=zeros(Tree.nNodes,1); chainok=1;
for ss=1:Tree.nSeg
    sn=Tree.segnodes{ss};
    covered(sn)=covered(sn)+1;
    chainok=chainok*all(Tree.nodenchild(sn(1:end-1))==1)*all(Tree.nodeparent(sn(2:end))==sn(1:end-1)')*all(diff(Tree.nodedate(sn))==1);
    if Tree.segparent(ss)>0
        chainok=chainok*(Tree.nodeseg(Tree.segparent(ss))<ss)*(Tree.nodenchild(Tree.segparent(ss))>1);
    end
end
fprintf('Tree (c): every node is in exactly one segment, should be one: %i \n',all(covered==1))
fprintf('Tree (c): each segment is a chain, starts after a branching node, and comes after its parent segment, should be one: %i \n',chainok)

%% (d) 3-D pi_S: identical slices give the same tree as 2-D; an identity slice stops branching in that period
Tree2=UncertainTransPathTree(T,T1,T2,pi_S,1);
Tree3=UncertainTransPathTree(T,T1,T2,repmat(pi_S,1,1,T2-T1+1),1);
fprintf('Tree (d): a 3-D pi_S with identical slices gives the same tree as the 2-D one, should be one: %i \n',isequal(Tree2.histories,Tree3.histories) && isequal(Tree2.nodeprob,Tree3.nodeprob) && isequal(Tree2.nodeS,Tree3.nodeS))
pi_S3=repmat(pi_S,1,1,T2-T1+1);
pi_S3(:,:,2)=eye(2);
Tree3=UncertainTransPathTree(T,T1,T2,pi_S3,1);
fprintf('Tree (d): with the middle slice the identity there are 2^2 histories, should be one: %i \n',Tree3.nHist==4)
fprintf('Tree (d): and S does not change in that period, should be one: %i \n',all(Tree3.historyS(:,T1+1)==Tree3.historyS(:,T1)))

%% (e) zero probabilities prune: pi_S=eye(2) is a single history
Tree=UncertainTransPathTree(T,T1,T2,eye(2),2);
fprintf('Tree (e): pi_S=eye(2) gives one history, at S_initial throughout, should be one: %i \n',Tree.nHist==1 && all(Tree.historyS==2))

%% (f) T1=1: S_1 is drawn, the period 1 nodes have no parent
Tree=UncertainTransPathTree(T,1,2,[0.5,0.5;0.5,0.5],1);
fprintf('Tree (f): with T1=1 there are two period 1 nodes, both without a parent, should be one: %i \n',sum(Tree.nodedate==1)==2 && all(Tree.nodeparent(Tree.nodedate==1)==0))
fprintf('Tree (f): and two segments without a parent segment, should be one: %i \n',sum(Tree.segparent==0)==2)
fprintf('Tree (f): and 4 histories, should be one: %i \n',Tree.nHist==4)

%% (g) prices and parameters on the tree, in each form they can be given
Tree=UncertainTransPathTree(T,T1,T2,pi_S,1);
xS=[1:T; 1000+(1:T)];     % n_S-by-T, value 1000*(S-1)+t
x1=1:T;                   % 1-by-T
xN=1:Tree.nNodes;         % 1-by-nNodes
PricePathG.a=xS; PricePathG.b=xS'; PricePathG.c=x1; PricePathG.d=x1'; PricePathG.e=xN;
ParamPathG.f=xS; ParamPathG.g=x1;
[~,PricePathNodes,ParamPathNodes]=UncertainTransPathTree(T,T1,T2,pi_S,1,PricePathG,ParamPathG);
xhand=(1000*(Tree.nodeS-1)+Tree.nodedate)';
fprintf('Tree (g): an n_S-by-T price is read at (S,t) of each node, either orientation, this should be zero: %.3e \n',max(abs([PricePathNodes.a-xhand,PricePathNodes.b-xhand])))
fprintf('Tree (g): a 1-by-T price is read at t, either orientation, this should be zero: %.3e \n',max(abs([PricePathNodes.c-Tree.nodedate',PricePathNodes.d-Tree.nodedate'])))
fprintf('Tree (g): a 1-by-nNodes price is one per node, this should be zero: %.3e \n',max(abs(PricePathNodes.e-xN)))
fprintf('Tree (g): ParamPath likewise, n_S-by-T and 1-by-T, this should be zero: %.3e \n',max(abs([ParamPathNodes.f-xhand,ParamPathNodes.g-Tree.nodedate'])))

output=struct();

end
