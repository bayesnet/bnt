function marginal = marginal_nodes(engine, nodes, t, fam)
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

if t < engine.T
  slice = t+1;
  nodes2 = nodes;
else % earliest t is T, so all nodes fit in one slice
  slice = engine.T;
  nodes2 = nodes + ss;
end
  
c = clq_containing_nodes(engine.jtree_engine, nodes2, fam);
assert(c >= 1);

%disp(['computing marginal on ' num2str(nodes) ' t = ' num2str(t)]);
%disp(['using ' num2str(nodes2) ' slice = ' num2str(slice) 'clq = ' num2str(c)]);

bigpot = engine.clpot{c, slice};

pot = marginalize_pot(bigpot, nodes2, engine.maximize);
marginal = pot_to_marginal(pot);

% we convert the domain to the unrolled numbering system
% so that update_ess extracts the right evidence.
marginal.domain = nodes+(t-1)*ss;
