function marginal = marginal_nodes(engine, query, add_ev)
% MARGINAL_NODES Compute the marginal on the specified query nodes (stab_cond_gauss)
% marginal = marginal_nodes(engine, query, add_ev)
%
% 'query' must be a singleton set.
% add_ev is an optional argument; if 1, we will "inflate" the marginal of observed nodes
% to their original size, adding 0s to the positions which contradict the evidence

if nargin < 3, add_ev = 0; end
if isempty(engine.evidence)
    hquery = query;  
else   
    hquery = [];
    for i = query
        if isempty(engine.evidence{i})
        hquery = [hquery i];
        end
    end
end

bnet = bnet_from_engine(engine);

nclq = length(engine.cliques);
clique = 0;
for i = 1:nclq
  if mysubset(hquery, engine.cliques{i})
    pot = struct(engine.clpot{i});
    %if mysubset(hquery, pot.cheaddom) | mysubset(hquery, pot.ddom)
    if mysubset(hquery, pot.domain)
     clique = i;
      break;
    end
  end
end

if isempty(hquery)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % If all requested variables are observed, no query is necessary %
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     marginal.mu = [];
     marginal.Sigma = [];
     marginal.T = 1.0;
     marginal.domain = query;
else
    if clique == 0
        marginal = marginal_difclq_nodes(engine, hquery);
    else 
        marginal = marginal_singleclq_nodes(engine, clique, hquery);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Change the format of output, so that it is identical to the %
    % format obtained by the same request for the junction-tree   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    marginal.domain = query;
    bnet = bnet_from_engine(engine);
    dquery = myintersect(bnet.dnodes,hquery);
    ns = bnet.node_sizes(dquery);
    if length(ns) == 0
    marginal.T = 1;
    else
        if length(ns) == 1
            ns = [1 ns];
        end
        marginal.T = reshape(marginal.T,ns);
    end
end
if add_ev
  bnet = bnet_from_engine(engine);
  %marginal = add_ev_to_dmarginal(marginal, engine.evidence, bnet.node_sizes);
  marginal = add_evidence_to_gmarginal(marginal, engine.evidence, bnet.node_sizes, bnet.cnodes);
end






