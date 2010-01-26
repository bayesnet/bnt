function engine = ff_inf_engine(bnet)
% FF_INF_ENGINE Factored frontier inference engine for DBNs
% engine = ff_inf_engine(bnet)
%
% The model must be topologically isomorphic to an HMM.
% In addition, each hidden node is assumed to have at most one observed child,
% and each observed child is assumed to have exactly one hidden parent.
%
% For details of this algorithm, see
%  "The Factored Frontier Algorithm for Approximate Inference in DBNs",
%   Kevin Murphy and Yair Weiss, UAI 2001.
%
% THIS IS HIGHLY EXPERIMENTAL CODE!

ss = length(bnet.intra);
onodes = bnet.observed;
hnodes = mysetdiff(1:ss, onodes);

[persistent_nodes, transient_nodes] = partition_dbn_nodes(bnet.intra, bnet.inter);
assert(isequal(onodes, transient_nodes));
assert(isequal(hnodes, persistent_nodes));

engine.onodes = onodes;
engine.hnodes = hnodes;
engine.marginals = [];
engine.fwd = [];
engine.back = [];
engine.CPDpot = [];
engine.filter = [];

obschild = zeros(1,ss);
for i=engine.hnodes(:)'
  %ocs = myintersect(children(bnet.dag, i), onodes);
  ocs = children(bnet.intra, i);
  assert(length(ocs) <= 1);
  if length(ocs)==1
    obschild(i) = ocs(1);
  end
end  
engine.obschild = obschild;


engine = class(engine, 'ff_inf_engine', inf_engine(bnet));

