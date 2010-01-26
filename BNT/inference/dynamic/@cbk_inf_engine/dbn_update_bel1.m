function [engine, loglik] = dbn_update_bel1(engine, evidence)
% DBN_UPDATE_BEL1 Update  the initial belief state (bk)
% engine = dbn_update_bel1(engine, evidence)
%
% evidence{i} has the evidence on node i for slice 1

bnet = bnet_from_engine(engine);
ss = length(bnet.intra);
CPDpot = cell(1,ss);      
t = 1;
for n=1:ss
  fam = family(bnet.dag, n);
  e = bnet.equiv_class(n, 1);
  CPDpot{n} = convert_to_pot(bnet.CPD{e}, engine.pot_type, fam(:), evidence);
end

onodes = find(~isemptycell(evidence));

[clpot, loglik] = enter_soft_evidence(engine.sub_engine1, engine.clq_ass_to_node1, CPDpot, onodes, engine.pot_type);

C  = length(engine.clusters);
newbel = cell(1,C);
for c=1:C
  k = engine.clq_ass_to_cluster1(c);
  newbel{c} = marginalize_pot(clpot{k}, engine.clusters{c});
end

engine.bel = newbel;
engine.bel_clpot = clpot;
engine.slice1 = 1;
