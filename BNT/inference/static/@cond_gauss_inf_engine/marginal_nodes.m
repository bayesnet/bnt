function marginal = marginal_nodes(engine, query, add_ev)
% MARGINAL_NODES Compute the marginal on the specified query nodes (cond_gauss)
% marginal = marginal_nodes(engine, query, add_ev)
%
% 'query' must be a singleton set
% add_ev is an optional argument; if 1, we will "inflate" the marginal of observed nodes
% to their original size, adding 0s to the positions which contradict the evidence

if nargin < 3, add_ev = 0; end

if length(query) ~= 1
  error('cond_gauss_inf_engine can only handle marginal queries on single nodes')
end
j = query;
bnet = bnet_from_engine(engine);

if myismember(j, bnet.cnodes)
  if ~myismember(j, engine.onodes)
    [m, C] = collapse_mog(engine.mu{j}, engine.Sigma{j}, engine.T);    
    marginal.mu = m;
    marginal.Sigma = C;
    marginal.T = 1.0; % single mixture component
  else
    marginal.mu = engine.evidence{j};
    k = bnet.node_sizes(j);
    marginal.Sigma = zeros(k,k);
    marginal.T = 1.0; % since P(E|E)=1
  end
else
  marginal = pot_to_marginal(marginalize_pot(engine.joint_dmarginal, j));
  if add_ev
    marginal = add_ev_to_dmarginal(marginal, engine.evidence, bnet.node_sizes);
  end
end

marginal.domain = query;
