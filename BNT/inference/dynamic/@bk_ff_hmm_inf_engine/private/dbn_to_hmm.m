function [prior, transmat] = dbn_to_hmm(bnet)
% DBN_TO_HMM Compute the discrete HMM matrices from a simple DBN
% [prior, transmat] = dbn_to_hmm(bnet)

onodes = bnet.observed;
ss = length(bnet.intra);
evidence = cell(1,2*ss);
hnodes = mysetdiff(1:ss, onodes);
prior = multiply_CPTs(bnet, [], hnodes, evidence);
transmat = multiply_CPTs(bnet, hnodes, hnodes+ss, evidence);
%obsmat1 = multiply_CPTs(bnet, hnodes, onodes, evidence);
%obsmat = multiply_CPTs(bnet, hnodes+ss, onodes+ss, evidence);
%obsmat1 = obsmat if the observation matrices are tied across slices



%%%%%%%%%%%%

function mat = multiply_CPTs(bnet, pdom, cdom, evidence)

% MULTIPLY_CPTS Make a matrix Pr(Y|X), where X represents all the parents, and Y all the children
% We assume the children have no intra-connections.
%
% e.g., Consider the DBN with interconnectivity i->i', j->j',k', k->i',k'
% Then transition matrix = Pr(i,j,k -> i',j',k') = Pr(i,k->i') Pr(j->j') Pr(j,k->k')

dom = [pdom cdom];
ns = bnet.node_sizes;
bigpot = dpot(dom, ns(dom));
for j=cdom(:)'
  e = bnet.equiv_class(j);
  fam = family(bnet.dag, j);
  pot = convert_to_pot(bnet.CPD{e}, 'd', fam(:), evidence);
  bigpot = multiply_by_pot(bigpot, pot);
end
psize = prod(ns(pdom));
csize = prod(ns(cdom));
T = pot_to_marginal(bigpot);
mat = reshape(T.T, [psize csize]);

          
