function engine = distribute_evidence(engine, root)

if isempty(engine.preorder{root})
  % this is the first time we have distributed from this root
  % memoize the order
  [jtree, preorder, postorder] = mk_rooted_tree(engine.jtree, root);
  preorder_children = cell(1,length(preorder));
  for n=preorder
    preorder_children{n} = children(jtree, n);
  end
  engine.preorder{root} = preorder;
  engine.preorder_children{root} = preorder_children;
else
  preorder = engine.preorder{root};
  preorder_children = engine.preorder_children{root};
end


% distribute from root (node to children)
for n=preorder(:)'
  for c=preorder_children{n}(:)'
    engine.clpot{c} = divide_by_pot(engine.clpot{c}, engine.seppot{n,c}); 
    engine.seppot{n,c} = marginalize_pot(engine.clpot{n}, engine.separator{n,c}, engine.maximize);
    engine.clpot{c} = multiply_by_pot(engine.clpot{c}, engine.seppot{n,c});
  end
end
