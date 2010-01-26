function [pot, MEU] = marginal_nodes(engine, d)

C = length(cliques);
%clpot = init_clpot(limid, cliques, d, clq_ass_to_node);
clpot = init_clpot(limid, cliques, [], clq_ass_to_node);

% collect to root
if 1
  % HUGIN
  seppot = cell(C, C);    % separators are implicitely initialized to 1s
  for n=postorder{di}(1:end-1)
    for p=parents(rooted_jtree{di}, n)
      %clpot{p} = divide_by_pot(clpot{n}, seppot{p,n}); % dividing by 1 is redundant
      seppot{p,n} = marginalize_pot(clpot{n}, separator{p,n});
      clpot{p} = multiply_by_pot(clpot{p}, seppot{p,n});
    end
  end
else
  % Shafer-Shenoy
  msg = cell(C,C);
  for n=postorder{di}(1:end-1)
    for c=children(rooted_jtree{di}, n)
      clpot{n} = multiply_by_pot(clpot{n}, msg{c,n});
    end
    p = parents(rooted_jtree{di}, n);
    %msg{n,p} = marginalize_pot(clpot{n}, cliques{p});
    msg{n,p} = marginalize_pot(clpot{n}, separator{n,p});
  end
  root = clq_ass_to_node(d);
  n=postorder{di}(end);
  assert(n == root);
  for c=children(rooted_jtree{di}, n)
    clpot{n} = multiply_by_pot(clpot{n}, msg{c,n});
  end
end	

fam = family(limid.dag, d);
pot = marginalize_pot(clpot{root}, fam);

%%%%%%%
jpot = compute_joint_pot_limid(limid);
pot2 = marginalize_pot(jpot, fam);
assert(approxeq_pot(pot, pot2))
%%%%%%

[policy, score] = extract_policy(pot);

e = limid.equiv_class(d);
limid.CPD{e} = set_params(limid.CPD{e}, 'policy', policy);

  
    
