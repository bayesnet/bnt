function engine = frontier_inf_engine(bnet)
% FRONTIER_INF_ENGINE Inference engine for DBNs which which uses the frontier algorithm.
% engine = frontier_inf_engine(bnet)
%
% The frontier algorithm extends the forwards-backwards algorithm to DBNs in the obvious way,
% maintaining a joint distribution (frontier) over all the nodes in a time slice.
% When all the hidden nodes in the DBN are persistent (have children in the next time slice),
% its theoretical running time is often similar to that of the junction tree algorithm,
% although in practice, this algorithm seems to very slow (at least in matlab).
% However, it is extremely simple to describe and implement.
%
% Suppose there are n binary nodes per slice, so the frontier takes O(2^n) space.
% Each time step takes between O(n 2^{n+1}) and O(n 2^{2n}) operations, depending on the graph structure.
% The lower bound is achieved by a set of n independent chains, as in a factorial HMM.
% The upper bound is achieved by a set of n fully interconnected chains, as in an HMM.
%
% The factor of n arises because we need to multiply in each CPD from slice t+1.
% The second factor depends on the size of the frontier to which we add the new node.
% In an FHMM, once we have added X(i,t+1), we can marginalize out X(i,t) from the frontier, since
% no other nodes depend on it; hence the frontier never contains more than n+1 nodes.
% In a fully coupled HMM, we must leave X(i,t) in the frontier until all X(j,t+1) have been
% added; hence the frontier will contain 2*n nodes at its peak.
%
% For details, see
%   "The Factored Frontier Algorithm for Approximate Inference in DBNs",
%   Kevin Murphy and Yair Weiss, UAI 01.

ns = bnet.node_sizes_slice;
onodes = bnet.observed;
ns(onodes) = 1;
ss = length(bnet.intra);

[engine.ops, engine.fdom] = best_first_frontier_seq(ns, bnet.dag);
engine.ops1 = 1:ss;

engine.fwdback = [];
engine.fwd_frontier = [];
engine.back_frontier = [];

engine.fdom1 = cell(1,ss);
for s=1:ss
  engine.fdom1{s} = 1:s;
end

engine = class(engine, 'frontier_inf_engine', inf_engine(bnet));


%%%%%%%%%

function [ops, frontier_set] = best_first_frontier_seq(ns, dag)
% BEST_FIRST_FRONTIER_SEQ Do a greedy search for the sequence of additions/removals to the frontier.
% [ops, frontier_set] = best_first_frontier_seq(ns, dag)
%
% We maintain 3 sets: the frontier (F), the right set (R), and the left set (L).
% The invariant is that the nodes in R are d-separated from L given F.
% We start with slice 1 in F and slice 2 in R.
% The goal is to move slice 1 from F to L, and slice 2 from R to F, so as to minimize the size
% of the frontier at each step, where the size(F) = product of the node-sizes of nodes in F.
% A node may be removed (from F to L) if it has no children in R.
% A node may be added (from R to F) if its parents are in F.
%
% ns(i) = num. discrete values node i can take on (i=1..ss, where ss = slice size)
% dag is the (2*ss) x (2*ss) adjacency matrix for the 2-slice DBN.

% Example:
%
% 4    9
% ^    ^
% |    |
% 2 -> 7
% ^    ^
% |    |
% 1 -> 6
% |    |
% v    v
% 3 -> 8
% |    |
% v    V
% 5    10
%
% ops = -4, -5, 6, -1, 7, -2, 8, -3, 9, 10

ss = length(ns);
ns = [ns(:)' ns(:)'];
ops = zeros(1,ss);
L = []; F = 1:ss; R = (1:ss)+ss;
frontier_set = cell(1,2*ss);
for s=1:2*ss
  remcost = inf*ones(1,2*ss);
  %disp(['L: ' num2str(L) ', F: ' num2str(F) ', R: ' num2str(R)]);
  maybe_removable = myintersect(F, 1:ss);
  for n=maybe_removable(:)'
    cs = children(dag, n);
    if isempty(myintersect(cs, R))
      remcost(n) = prod(ns(mysetdiff(F, n)));
    end
  end
  %remcost
  if any(remcost < inf)
    n = argmin(remcost);
    ops(s) = -n;
    L = myunion(L, n);
    F = mysetdiff(F, n);
  else
    addcost = inf*ones(1,2*ss);
    for n=R(:)'
      ps = parents(dag, n);
      if mysubset(ps, F)
	addcost(n) = prod(ns(myunion(F, [ps n])));
      end
    end
    %addcost
    assert(any(addcost < inf));
    n = argmin(addcost);
    ops(s) = n;
    R  = mysetdiff(R, n);
    F = myunion(F, n);
  end
  %fprintf('op at step %d = %d\n\n', s, ops(s));
  frontier_set{s} = F;
end
