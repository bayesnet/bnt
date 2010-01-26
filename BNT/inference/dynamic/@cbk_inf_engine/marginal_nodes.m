function marginal = marginal_nodes(engine, nodes, t, fam)
% MARGINAL_NODES Compute the marginal on the specified query nodes (bk)
%
%   marginal = marginal_nodes(engine, i, t)
% returns Pr(X(i,t) | Y(1:T)), where X(i,t) is the i'th node in the t'th slice.
% If enter_evidence used filtering instead of smoothing, this will return  Pr(X(i,t) | Y(1:t)).
%
%   marginal = marginal_nodes(engine, query, t)
% returns Pr(X(query(1),t), ... X(query(end),t) | Y(1:T)),
% where X(q,t) is the q'th node in the t'th slice. If q > ss (slice size), this is equal
% to X(q mod ss, t+1). That is, 't' specifies the time slice of the earliest node.
% 'query' cannot span more than 2 time slices.
% Example:
% Consider a DBN with 2 nodes per slice.
% Then t=2, nodes=[1 3] refers to node 1 in slice 2 and node 1 in slice 3.

if nargin < 3, t = 1; end
if nargin < 4, fam = 0; else fam = 1; end


% clpot{t} contains slice t-1 and t
% Example
% clpot #: 1    2    3
% slices:  1  1,2  2,3
% For filtering, we must take care not to take future evidence into account.
% For smoothing, clpot{1} does not exist.

bnet = bnet_from_engine(engine);
ss = length(bnet.intra);

nodes2 = nodes;
if ~engine.filter
  if t < engine.T
    slice = t+1;
  else % earliest t is T, so all nodes fit in one slice
    slice = engine.T;
    nodes2 = nodes + ss;
  end
else
  if t == 1
   slice = 1;
  else
    if all(nodes<=ss)
      slice = t;
      nodes2 = nodes + ss;
    elseif t == engine.T
      slice = t;
    else
      slice = t + 1;
    end
  end
end
  
if engine.filter & t==1
  c = clq_containing_nodes(engine.sub_engine1, nodes2, fam);
else
  c = clq_containing_nodes(engine.sub_engine, nodes2, fam);
end
assert(c >= 1);
bigpot = engine.clpot{c, slice};

pot = marginalize_pot(bigpot, nodes2);
marginal = pot_to_marginal(pot);

% we convert the domain to the unrolled numbering system
% so that update_ess extracts the right evidence.
marginal.domain = nodes+(t-1)*ss;
