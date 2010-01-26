function [engine, loglik] = initialize_engine(engine)
%initialize
bnet = bnet_from_engine(engine);
ns = bnet.node_sizes(:);
N = length(bnet.dag);

pot_type = 'scg'
check_for_cd_arcs([], bnet.cnodes, bnet.dag);

% Evaluate CPDs with evidence, and convert to potentials  
pot = cell(1, N);
C = length(engine.cliques);
inited = zeros(1, C);
clpot = cell(1, C);
evidence = cell(1, N);
for n=1:N
  fam = family(bnet.dag, n);
  e = bnet.equiv_class(n);
  pot{n} = CPD_to_scgpot(bnet.CPD{e}, fam, ns, bnet.cnodes, evidence);
  cindex = engine.clq_ass_to_node(n);
  if inited(cindex)
      %clpot{cindex} = direct_combine_pots(clpot{cindex}, pot{n});
      clpot{cindex} = direct_combine_pots(pot{n}, clpot{cindex});
  else
      clpot{cindex} = pot{n};
      inited(cindex) = 1;
  end
end

for i=1:C
    if inited(i) == 0
        clpot{i} = scgpot([], [], [], []);
    end
end

seppot = cell(C, C);
% separators are is not need to initialize

% collect to root (node to parents)
for n=engine.postorder(1:end-1)
  for p=parents(engine.jtree, n)
      [margpot, comppot] = complement_pot(clpot{n}, engine.separator{p,n});
      margpot = marginalize_pot(clpot{n}, engine.separator{p,n});
      clpot{n} = comppot;
      %seppot{p, n} = margpot;
      clpot{p} = combine_pots(clpot{p}, margpot);
      %clpot{p} = combine_pots(margpot, clpot{p});
  end
end

temppot = clpot;
%temppot = clpot{engine.root};
for n=engine.preorder
  for c=children(engine.jtree, n)
    seppot{n,c} = marginalize_pot(temppot{n}, engine.separator{n,c});
    %seppot{n,c} = marginalize_pot(clpot{n}, engine.separator{n,c});
    %clpot{c} = direct_combine_pots(clpot{c}, seppot{n,c});
    temppot{c} = direct_combine_pots(temppot{c}, seppot{n,c});
  end
end

engine.clpot = clpot;
engine.seppot = seppot;


