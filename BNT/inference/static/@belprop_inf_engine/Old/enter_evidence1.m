function engine = enter_evidence(engine, evidence)

doms = engine.fgraph.doms;
ndoms = length(doms);
ns = engine.fgraph.node_sizes;
obs = find(~isemptycell(evidence));
cobs = myintersect(obs, engine.fgraph.cnodes);
dobs = myintersect(obs, engine.fgraph.dnodes);
ns(cobs) = 0;
ns(dobs) = 1;

% recompute the weight of each domain now that we know what nodes are observed
for i=1:ndoms
  engine.dom_weight(i) = prod(ns(engine.fgraph.doms{i}));
end

% prime each local kernel with evidence (if any)
local_kernel = cell(1, ndoms);
for i=1:length(engine.fgraph.kernels_of_type)
  u = engine.fgraph.kernels_of_type{i};
  local_kernel(u) = kernel_to_dpots(engine.fgraph.kernels{i}, evidence, engine.fgraph.domains_of_type{i});
end
  
% initialise all msgs to 1s
msg = cell(ndoms, ndoms);
for i=1:ndoms
  nbrs = engine.fgraph.nbrs{i};
  for j=nbrs(:)'
    dom = engine.fgraph.sepset{i,j};
    msg{i,j} = dpot(dom, ns(dom));
  end
end

prod_of_msg = cell(1, ndoms);
bel = cell(1, ndoms);
old_bel = cell(1, ndoms);

converged = 0;
iter = 1;
while ~converged & (iter <= engine.max_iter)
  
  % each node multiplies all its incoming msgs
  for i=1:ndoms
    prod_of_msg{i} = dpot(doms{i}, ns(doms{i}));
    nbrs = engine.fgraph.nbrs{i};
    for j=nbrs(:)'
      prod_of_msg{i} = multiply_by_pot(prod_of_msg{i}, msg{j,i});
    end
  end
  
  % each node computes its local belief
  old_bel = bel;
  for i=1:ndoms
    bel{i} = normalize_pot(multiply_pots(prod_of_msg{i}, local_kernel{i}));
  end

  % converged?
  if iter==1
    converged = 0;
  else
    converged = 1;
    for i=1:ndoms
      belT = get_params(bel{i}, 'table');
      old_belT = get_params(old_bel{i}, 'table');
      if ~approxeq(belT, old_belT, engine.tol)
	converged = 0;
	break;
      end
    end
  end

  if ~converged
    old_msg = msg;
    % each node sends a msg to each of its neighbors
    for i=1:ndoms
      nbrs = engine.fgraph.nbrs{i};
      for j=nbrs(:)'
	% multiply all incoming msgs except from j
	temp = prod_of_msg{i};
	temp = divide_by_pot(temp, old_msg{j,i});
	% send msg from i to j
	temp = multiply_by_pot(temp, local_kernel{i});
	msg{i,j} = normalize_pot(marginalize_pot(temp, engine.fgraph.sepset{i,j}));
      end
    end
  end

  iter = iter + 1
end

engine.marginal_domains = bel;
%for i=1:ndoms  
  %engine.marginal_domains{i} = get_params(bel{i}, 'table');
%end
