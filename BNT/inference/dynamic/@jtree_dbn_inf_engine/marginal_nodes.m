function marginal = marginal_nodes(engine, nodes, t, add_ev, fam)
% MARGINAL_NODES Compute the marginal on the specified query nodes (bk)
%
%   marginal = marginal_nodes(engine, i, t)
% returns Pr(X(i,t) | Y(1:T)), where X(i,t) is the i'th node in the t'th slice.
%
%   marginal = marginal_nodes(engine, query, t)
% returns Pr(X(query(1),t), ... X(query(end),t) | Y(1:T)),
% where 't' specifies the time slice of the earliest node in the query.
% 'query' cannot span more than 2 time slices.
%
% Example:
% Consider a DBN with 2 nodes per slice.
% Then t=2, nodes=[1 3] refers to node 1 in slice 2 and node 1 in slice 3.
%
% marginal = marginal_nodes(engine, nodes, t, add_ev, fam)
% add_ev is an optional argument; if 1, we will "inflate" the marginal of observed nodes
% to their original size, adding 0s to the positions which contradict the evidence
   
if nargin < 3, t = 1; end
if nargin < 4, add_ev = 0; end
if nargin < 5, fam = 0; end

bnet = bnet_from_engine(engine);
ss = length(bnet.intra);

if t==1 | t==engine.T
  slice = t;
  nodes2 = nodes;
elseif mysubset(nodes, engine.persist)
  slice = t-1;
  nodes2 = nodes+ss;
else
  slice = t;
  nodes2 = nodes;
end

%disp(['computing marginal on ' num2str(nodes) ' t = ' num2str(t) ' fam = ' num2str(fam)]);

if t==engine.T
  c = clq_containing_nodes(engine.jtree_engine1, nodes2, fam);
else
  c = clq_containing_nodes(engine.jtree_engine, nodes2, fam);
end
if c == -1
  error(['no clique contains ' nodes2])
end


%disp(['using ' num2str(nodes2) ' slice = ' num2str(slice) ' clq = ' num2str(c)]);

bigpot = engine.clpot{c, slice};

pot = marginalize_pot(bigpot, nodes2, engine.maximize);
%pot = normalize_pot(pot);
marginal = pot_to_marginal(pot);


% we convert the domain to the unrolled numbering system
% so that add_ev_to_dmarginal (maybe called in update_ess) extracts the right evidence.
marginal.domain = nodes+(t-1)*ss;

if add_ev
  marginal = add_ev_to_dmarginal(marginal, engine.evidence, engine.node_sizes);
end    

