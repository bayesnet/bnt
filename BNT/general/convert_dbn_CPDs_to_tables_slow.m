function CPDpot = convert_dbn_CPDs_to_tables_slow(bnet, evidence)
% CONVERT_DBN_CPDS_TO_TABLES_SLOW Convert CPDs of (possibly instantiated) DBN nodes to tables
% CPDpot = convert_dbn_CPDs_to_tables_slow(bnet, evidence)
%
% CPDpot{n,t} is a table containing P(n,t|pa(n,t), ev)
% All hidden nodes are assumed to be discrete
%
% Non-vectorized method; this is less efficient for long sequences of observed Gaussian
% nodes, because of the (unnecessary) repeated matrix inversion.

obs_bitv = ~isemptycell(evidence(:));
[ss T] = size(evidence);
ns = bnet.node_sizes(:);

CPDpot = cell(ss,T); 

t = 1;
for n=1:ss
  %ps = engine.bnet_parents{n};
  ps = parents(bnet.dag, n);
  e = bnet.equiv_class(n, 1);
  if ~any(obs_bitv(ps))
    CPDpot{n,t} = convert_CPD_to_table_hidden_ps(bnet.CPD{e}, evidence{n,t});
  else
    CPDpot{n,t} = convert_to_table(bnet.CPD{e}, [ps n], evidence(:,1));
  end
end
for t=2:T
  for n=1:ss
    self = n+ss;
    ps = parents(bnet.dag, self);
    e = bnet.equiv_class(n, 2);
    if ~any(obs_bitv(ps))
      CPDpot{n,t} = convert_CPD_to_table_hidden_ps(bnet.CPD{e}, evidence{n,t});
    else
      CPDpot{n,t} = convert_to_table(bnet.CPD{e}, [ps self], evidence(:,t-1:t));
    end
  end
end       


