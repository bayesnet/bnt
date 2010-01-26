function [bel, niter] = parallel_protocol(engine, evidence, pot_type, local_kernel, msg)

bnet = bnet_from_engine(engine);
ns = bnet.node_sizes;
onodes = find(~isemptycell(evidence));

ndoms = length(engine.gdl.doms);
prod_of_msg = cell(1, ndoms);
bel = cell(1, ndoms);
old_bel = cell(1, ndoms);

converged = 0;
iter = 1;
while ~converged & (iter <= engine.max_iter)
  
  % each node multiplies all its incoming msgs and computes its local belief
  old_bel = bel;
  for i=1:ndoms
    prod_of_msg{i} = mk_initial_pot(pot_type, engine.gdl.doms{i}, ns, bnet.cnodes, onodes);
    nbrs = engine.gdl.nbrs{i};
    for j=nbrs(:)'
      prod_of_msg{i} = multiply_by_pot(prod_of_msg{i}, msg{j,i});
    end
    bel{i} = normalize_pot(multiply_by_pot(local_kernel{i}, prod_of_msg{i}));
  end

  if ~isempty(engine.fid)
    for i=1:ndoms
      tmp = pot_to_marginal(bel{i});
      %fprintf(engine.fid, '%9.7f ', tmp.T(1));
      fprintf(engine.fid, '%9.7f ', tmp.U(1));
    end
    %fprintf(engine.fid, '  U ');
    %for i=1:ndoms
    %  tmp = pot_to_marginal(bel{i});
    %  fprintf(engine.fid, '%9.7f ', tmp.U(1));
    %end
    fprintf(engine.fid, '\n');
  end

  % converged?
  if iter==1
    converged = 0;
  else
    converged = 1;
    for i=1:ndoms
      if ~approxeq_pot(bel{i}, old_bel{i}, engine.tol)
	converged = 0;
	break;
      end
    end
  end

  if ~converged
    old_msg = msg;
    % each node sends a msg to each of its neighbors
    for i=1:ndoms
      nbrs = engine.gdl.nbrs{i};
      for j=nbrs(:)'
	% multiply all incoming msgs except from j
	temp = prod_of_msg{i};
	temp = divide_by_pot(temp, old_msg{j,i});
	% send msg from i to j
	temp = multiply_by_pot(temp, local_kernel{i});
	temp2 = marginalize_pot(temp, engine.gdl.sepset{i,j}, engine.maximize);
	msg{i,j} = normalize_pot(temp2);
      end
    end
  end

  iter = iter + 1;
end


niter = iter-1;

if 0
for i=1:ndoms
  prod_of_msg{i} = mk_initial_pot(pot_type, engine.gdl.doms{i}, ns, bnet.cnodes, onodes);
  nbrs = engine.gdl.nbrs{i};
  for j=nbrs(:)'
    prod_of_msg{i} = multiply_by_pot(prod_of_msg{i}, msg{j,i});
  end
  bel{i} = normalize_pot(multiply_by_pot(local_kernel{i}, prod_of_msg{i}));
end
end
