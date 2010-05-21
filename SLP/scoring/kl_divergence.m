function KLdiv = KL_divergence(bnetP, bnetQ)
% KL_DIVERGENCE computes the Kullback-Leibler divergence between two BNET distributions
% KLdiv = KL_divergence(bnetP, bnetQ)
%
% Output :
%   div = sum_x  P(x).log(P(x)/Q(x))
%
% Rem : 
%   This version is optimized for speed, but can use too many memory
%     ==> if you have a memory problem, use kl_divergence2 instead
%

%   ONLY FOR TABULAR NODES
%   Make sure that you have done the params learning.
%  
%   V1.1 : 8 oct 2004 (Ph. Leray - philippe.leray@univ-nantes.fr)

N = size(bnetP.dag,1);
N2 = size(bnetQ.dag,1);
ns= bnetP.node_sizes;
ns2= bnetQ.node_sizes;
if N~=N2, error('size of dags must be the same'), end
if ns~=ns2, error('node sizes of dags must be the same'), end
tiny = exp(-700);
KLdiv=0;

inst = ind2subv(ns, 1:prod(ns)); 
  %Px=1; Qx=1;
  for i=1:N,
    ps = parents(bnetP.dag, i);
    %e = bnetP.equiv_class(i);
    %[tmp Px(:,i)] = prob_node(bnetP.CPD{e}, inst(:,i)', inst(:,ps)');
    [tmp Px(:,i)] = prob_node(bnetP.CPD{i}, inst(:,i)', inst(:,ps)');

    ps = parents(bnetQ.dag, i);
    %e = bnetQ.equiv_class(i);
    %[tmp Qx(:,i)] = prob_node(bnetQ.CPD{e}, inst(:,i)', inst(:,ps)');
    [tmp Qx(:,i)] = prob_node(bnetQ.CPD{i}, inst(:,i)', inst(:,ps)');
  end
 Px=prod(Px,2);
 Px = Px + (Px==0)*tiny; % replace 0s by tiny
 Qx=prod(Qx,2);
 Qx = Qx + (Qx==0)*tiny; % replace 0s by tiny

 %%%%% Faut-il diviser par le nb de configurations possibles ? (sum => mean)
 KLdiv = sum(Px.*log(Px./Qx));
  %end

