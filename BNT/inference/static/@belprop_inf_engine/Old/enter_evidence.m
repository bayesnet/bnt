function engine = enter_evidence(engine, evidence)

doms = engine.fg.doms;
ndoms = length(doms);
ns = engine.fg.node_sizes;
obs = find(~isemptycell(evidence));
cobs = myintersect(obs, engine.fg.cnodes);
dobs = myintersect(obs, engine.fg.dnodes);
ns(cobs) = 0;
ns(dobs) = 1;

% prime each local kernel with evidence (if any)
local_kernel = cell(1, ndoms);
for i=1:length(engine.fg.kernels_of_type)
  u = engine.fg.kernels_of_type{i};
  local_kernel(u) = kernel_to_dpots(engine.fg.kernels{i}, evidence, engine.fg.domains_of_type{i});
end
  
% initialise all msgs to 1s
nedges = engine.fg.nedges;
msg = cell(1, nedges);
for i=1:nedges
  msg{i} = dpot(engine.fg.sepset{i}, ns(engine.fg.sepset{i}));
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
    nbrs = engine.fg.nbrs{i};
    for j=1:length(nbrs)
      ndx = engine.fg.edge_ndx(j,i);
      prod_of_msg{i} = multiply_by_pot(prod_of_msg{i}, msg{ndx});
    end
  end
  old_msg = msg;
  
  % each node computes its local belief
  for i=1:ndoms
    bel{i} = normalize_pot(multiply_pots(prod_of_msg{i}, local_kernel{i}));
  end

  % converged?
  converged = 1;
  for i=1:ndoms
    if ~approxeq(bel{i}, old_bel{i}, engine.tol)
      converged = 0;
      break;
    end
  end

  if ~converged
    % each node sends a msg to each of its neighbors
    for i=1:ndoms
      nbrs = engine.fg.nbrs{i};
      for j=1:length(nbrs)
	% multiply all incoming msgs except from j
	temp = prod_of_msg{i};
	ndx = engine.fg.edge_ndx(j,i);
	temp = divide_by_pot(temp, old_msg{ndx});
	% send msg from i to j
	temp = multiply_by_pot(temp, local_kernel{i});
	ndx = engine.fg.edge_ndx(i,j);
	msg{ndx} = normalize_pot(marginalize_pot(temp, engine.fg.sepset{ndx}));
      end
    end
  end

  iter = iter + 1;
end

  
engine.marginal = bel;
