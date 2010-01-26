function marginal = marginal_nodes(engine, query, add_ev)
% MARGINAL_NODES Compute the marginal on the specified query nodes (loopy)
% marginal = marginal_nodes(engine, query, add_ev)
%
% 'query' must be a single node.
% add_ev is an optional argument; if 1, observed nodes will be set to their original size,
% otherwise they will be treated like points.
   
if nargin < 3, add_ev = 0; end

if length(query) > 1
  error('can only compute marginal on single nodes or families')
end
bnet = bnet_from_engine(engine);
ns = bnet.node_sizes(:);

switch engine.msg_type
 case 'd',
  T = engine.marginal{query};
  if ~add_ev
    marginal.T = shrink_obs_dims_in_table(T, query, engine.evidence);
  else
    marginal.T = T;
  end
  marginal.domain = query;
 
 case 'g',
  if engine.disconnected_nodes_bitv(query)
    marginal.T = 1;
    marginal.domain = query;
    if add_ev
      marginal = add_ev_to_dmarginal(marginal, engine.evidence, ns)
    end
    return;
  end

  marginal = engine.marginal{query};
  marginal.domain = query;
  if ~add_ev
    marginal = shrink_obs_dims_in_gaussian(marginal, query, engine.evidence, ns);
  end
end

