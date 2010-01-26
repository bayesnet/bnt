function engine = collect_evidence(engine, root)

if isempty(engine.postorder{root})
  % this is the first time we have collected to this root
  % memoize the order
  [jtree, preorder, postorder] = mk_rooted_tree(engine.jtree, root);
  postorder_parents = cell(1,length(postorder));
  for n=postorder(1:end-1)
    postorder_parents{n} = parents(jtree, n);
  end
  engine.postorder{root} = postorder;
  engine.postorder_parents{root} = postorder_parents;
else
  postorder = engine.postorder{root};
  postorder_parents = engine.postorder_parents{root};
end

C = length(engine.clpot);
seppot = cell(C, C);
% separators are implicitely initialized to 1s

% collect to root (node to parents)
for n=postorder(1:end-1)
  for p=postorder_parents{n}
    %clpot{p} = divide_by_pot(clpot{n}, seppot{p,n}); % dividing by 1 is redundant
    engine.seppot{p,n} = marginalize_pot(engine.clpot{n}, engine.separator{p,n}, engine.maximize);
    engine.clpot{p} = multiply_by_pot(engine.clpot{p}, engine.seppot{p,n});
  end
end
