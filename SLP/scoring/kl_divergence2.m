function KLdiv = KL_divergence2(bnetP, bnetQ)
% KL_DIVERGENCE2 computes the Kullback-Leibler divergence between two BNET distributions
% KLdiv = KL_divergence2(bnetP, bnetQ)
%
% Output :
%   div = sum_x  P(x).log(P(x)/Q(x))
%
% Rem : 
%   This version is optimized for memory use, but quite slow !!!
%     ==> if you have no memory problem, use kl_divergence instead
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

for i=1:prod(ns),
  inst = ind2subv(ns, i); % i'th instantiation
  Px=1; Qx=1;
  for i=1:N,
    ps = parents(bnetP.dag, i);
    e = bnetP.equiv_class(i);
    [tmp Pxi] = prob_node(bnetP.CPD{e}, inst(i), inst(ps)');
    Px=Px*Pxi;
    ps = parents(bnetQ.dag, i);
    e = bnetQ.equiv_class(i);
    [tmp Qxi] = prob_node(bnetQ.CPD{e}, inst(i), inst(ps)');
    Qx=Qx*Qxi;
  end
    Px = Px + (Px==0)*tiny; % replace 0s by tiny
    Qx = Qx + (Qx==0)*tiny; % replace 0s by tiny
  KLdiv = KLdiv + Px*log(Px/Qx);
end

