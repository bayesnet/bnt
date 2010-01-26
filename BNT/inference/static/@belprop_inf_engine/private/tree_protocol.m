function bel = tree_protocol(engine, evidence, pot_type, local_kernel, msg)

bnet = bnet_from_engine(engine);
ns = bnet.node_sizes;
onodes = find(~isemptycell(evidence));

ndoms = length(engine.gdl.doms);
prod_of_msg = cell(1, ndoms);
bel = cell(1, ndoms);
  
% collect to root (node to parents)
for n=engine.postorder
  % absorb msgs from children
  prod_of_msg{n} = mk_initial_pot(pot_type, engine.gdl.doms{n}, ns, bnet.cnodes, onodes);
  for c=children(engine.tree, n)
    prod_of_msg{n} = multiply_by_pot(prod_of_msg{n}, msg{c,n});
  end
  % send msg to parents
  for p=parents(engine.tree, n)
    if iter==1
      temp = prod_of_msg{n};
    else
      temp = divide_by_pot(prod_of_msg{n}, old_msg{p,n});
    end
    temp = multiply_by_pot(temp, local_kernel{n});
    temp2 = marginalize_pot(temp, engine.gdl.sepset{n,p}, engine.maximize);
    %fprintf('%d sends %d\n', n, p);
    msg{n,p} = normalize_pot(temp2);
  end
end

% distribute from root (node to children)
for n=engine.preorder
  % absorb from parents
  %prod_of_msg{n} = mk_initial_pot(pot_type, doms{n}, ns, cnodes, onodes);
  for p=parents(engine.tree, n)
    prod_of_msg{n} = multiply_by_pot(prod_of_msg{n}, msg{p,n});
  end
  bel{n} = normalize_pot(multiply_pots(prod_of_msg{n}, local_kernel{n}));
  % send msg to children
  for c=children(engine.tree, n)
    temp = divide_by_pot(prod_of_msg{n}, msg{c,n});
    temp = multiply_by_pot(temp, local_kernel{n});
    temp2 = marginalize_pot(temp, engine.gdl.sepset{n,c}, engine.maximize);
    %fprintf('%d sends %d\n', n, c);
    msg{n,c} = normalize_pot(temp2);
  end
end
