function [bnet, vals] = mk_minimal_qmr_bnet(G, inhibit, leak, prior, pos, neg, pos_only)
% MK_MINIMAL_QMR_BNET Make a QMR model which only contains the observed findings
% [bnet, vals] = mk_minimal_qmr_bnet(G, inhibit, prior, leak, pos, neg)
%
% Input:
% G(i,j) = 1 iff there is an arc from disease i to finding j
% inhibit(i,j) = inhibition probability on i->j arc
% leak(j) = inhibition prob. on leak->j arc
% prior(i) = prob. disease i is on
% pos = list of leaves that have positive observations
% neg = list of leaves that have negative observations
% pos_only = 1 means only include positively observed leaves in the model - the negative
%   ones are absorbed into the prior terms
%
% Output:
% bnet
% vals is their value

if pos_only
  obs = pos;
else
  obs = myunion(pos, neg);
end
Nfindings = length(obs);
[Ndiseases maxNfindings] = size(inhibit);
N = Ndiseases + Nfindings;
finding_node = Ndiseases+1:N;

% j = finding_node(i) means the i'th finding node is the j'th node in the bnet
% k = obs(i) means the i'th observed (positive) finding is the k'th finding overall
% If all findings are observed, and posonly = 0, we have i = obs(i) for all i.

%dag = sparse(N, N);
dag = zeros(N, N);
dag(1:Ndiseases, Ndiseases+1:N) = G(:,obs);

ns = 2*ones(1,N);
bnet = mk_bnet(dag, ns, 'observed', finding_node);

CPT = cell(1, Ndiseases);
for d=1:Ndiseases
  CPT{d} = [1-prior(d) prior(d)];
end

if pos_only
  % Fold in the negative evidence into the prior
  for i=1:length(neg)
    n = neg(i);
    ps = parents(G,n);
    for pi=1:length(ps)
      p = ps(pi);
      q = inhibit(p,n);
      CPT{p} = CPT{p} .* [1 q];
    end
    % Arbitrarily attach the leak term to the first parent
    p = ps(1);
    q = leak(n);
    CPT{p} = CPT{p} .* [q q];
  end
end

for d=1:Ndiseases
  bnet.CPD{d} = tabular_CPD(bnet, d, CPT{d}');
end

for i=1:Nfindings
  fnode = finding_node(i);
  fid = obs(i);
  ps = parents(G, fid);
  bnet.CPD{fnode} = noisyor_CPD(bnet, fnode, leak(fid), inhibit(ps, fid));
end

obs_nodes = finding_node;
vals = sparse(1, maxNfindings);
vals(pos) = 2;
vals(neg) = 1;
vals = full(vals(obs));





