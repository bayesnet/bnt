function pot_type = determine_pot_type(model, onodes, nodes)
% DETERMINE_POT_TYPE Determine the type of potential based on the evidence pattern.
% pot_type = determine_pot_type(model, onodes, nodes)
%
% If there are any utility nodes, pot_type = 'u'
% else
% If all hidden nodes are discrete, pot_type = 'd'.
% If all hidden nodes are continuous, pot_type = 'g' (Gaussian).
% If some hidden nodes are discrete, and some cts, pot_type = 'cg' (conditional Gaussian).
%
% nodes defaults to all nodes in graph

nnodes = length(model.node_sizes);
if nargin < 3, nodes = 1:nnodes; end

hnodes = mysetdiff(nodes, onodes);
if isfield(model, 'limid') %~isempty(model.utility_nodes)
  pot_type = 'u';
elseif isempty(myintersect(model.cnodes, hnodes))
  pot_type = 'd';
elseif mysubset(hnodes, model.cnodes)
  pot_type = 'g';
else
  pot_type = 'cg';
end
