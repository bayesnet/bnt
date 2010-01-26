function [engine, loglik] = dbn_update_bel(engine, evidence)
% DBN_UPDATE_BEL Update the belief state (bk)
% [engine, loglik] = dbn_update_bel(engine, evidence)
%
% evidence{i,1} contains the evidence on node i in slice t-1
% evidence{i,2} contains the evidence on node i in slice t

oldbel = engine.bel;

ss = size(evidence, 1);
bnet = bnet_from_engine(engine);
CPDpot = cell(1, ss);
for n=1:ss
  fam = family(bnet.dag, n, 2);
  e = bnet.equiv_class(n, 2);
  CPDpot{n} = convert_to_pot(bnet.CPD{e}, engine.pot_type, fam(:), evidence);
end

observed = ~isemptycell(evidence);
onodes2 = find(observed(:));
clqs = [engine.clq_ass_to_cluster(:,1); engine.clq_ass_to_node(:,2)];
pots = [oldbel(:); CPDpot(:)];

[clpot, loglik] = enter_soft_evidence(engine.sub_engine, clqs, pots, onodes2(:), engine.pot_type);

C = length(engine.clusters);
newbel = cell(1,C);
for c=1:C
  k = engine.clq_ass_to_cluster(c,2);
  cl = engine.clusters{c};
  newbel{c} = marginalize_pot(clpot{k}, cl+ss); % extract slice 2 posterior
  newbel{c} = set_domain_pot(newbel{c}, cl); % shift back to slice 1 for re-use as prior
end

engine.bel = newbel;
engine.bel_clpot = clpot;
engine.slice1 = 0;
