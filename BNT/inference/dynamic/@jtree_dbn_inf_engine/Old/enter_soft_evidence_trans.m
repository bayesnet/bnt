function [clpot, loglik] = enter_soft_evidence(engine, CPDpot, observed, pot_type)
% ENTER_SOFT_EVIDENCE Add the specified soft evidence to the network (jtree_dbn)
% [clpot, loglik] = enter_soft_evidence(engine, CPDpot, observed, pot_type, filter)

[ss T] = size(CPDpot);
Q = length(engine.jtree_struct.cliques);
clpot = cell(Q,T); % clpot{t} contains evidence from slices (t-1, t) 
seppot = cell(Q,Q,T);
ll = zeros(1,Q);
logscale = zeros(1,T);
bnet = bnet_from_engine(engine);

% Forwards pass.
% Compute distribution on clq C,
% where C is the out interface to (t-1,t).
% Then pass this to clq D, where D is the in inferface to (t+1,t).
% Then propagate from D to later slices.

C = engine.out_clq;
assert(C==engine.jtree_struct.root_clq);
D = engine.in_clq;
slice1 = 1:ss;
slice2 = slice1 + ss;
Ntransient = length(engine.transient);
trans = cell(Ntransient,1);
for t=1:T
  if t==1
    pots = [CPDpot(:,1); CPDpot(engine.persist, 2)];
    clqs = engine.jtree_struct.clq_ass_to_node([slice1 engine.persist+ss]);
    obs = find(observed(:,1:2));
  elseif t==T
    clqs = [D engine.jtree_struct.clq_ass_to_node(engine.transient)];
    phiC = set_domain_pot(phiC, engine.interface); % shift back to slice 1
    for i=1:Ntransient
      trans{i} = CPDpot{engine.transient(i), t};
      trans{i} = set_domain_pot(trans{i}, domain_pot(trans{i})-ss); % shift back to slice 1
    end
    pots = [ {phiC}; trans]; 
    obs = find(observed(:,T));
  else
    clqs = [D engine.jtree_struct.clq_ass_to_node([engine.transient engine.persist+ss])];
    phiC = set_domain_pot(phiC, engine.interface); % shift back to slice 1
    for i=1:Ntransient
      trans{i} = CPDpot{engine.transient(i), t};
      trans{i} = set_domain_pot(trans{i}, domain_pot(trans{i})-ss); % shift back to slice 1
    end
    pots = [ {phiC}; trans; CPDpot(engine.persist, t+1)]; 
    obs = find(observed(:,t:t+1));
  end
  [clpot(:,t), seppot(:,:,t)] =  init_pot(engine.jtree_struct.cliques, clqs, pots, pot_type, ...
					 obs, bnet.node_sizes(:), bnet.cnodes);
  [clpot(:,t), seppot(:,:,t)] = collect_evidence(clpot(:,t), seppot(:,:,t), engine.maximize, ...
						    engine.jtree_struct.postorder, ...
						    engine.jtree_struct.postorder_parents,...
						    engine.jtree_struct.separator);

  for c=1:Q
    [clpot{c,t}, ll(c)] = normalize_pot(clpot{c,t});
  end
  logscale(t) = ll(C);

  phiC = marginalize_pot(clpot{C,t}, engine.interface+ss, engine.maximize);
end



% Backwards pass.
% Pass evidence from clq C to clq D,
% where C is the in interface to (t,t+1) and D is the out inferface to (t-1,t)
% Then propagate evidence from D to earlier slices.
C = engine.in_clq;
D = engine.out_clq;
for t=T:-1:1
  [clpot(:,t), seppot(:,:,t)] = distribute_evidence(clpot(:,t), seppot(:,:,t), engine.maximize, ...
						    engine.jtree_struct.preorder, ...
						    engine.jtree_struct.preorder_children, ...
						    engine.jtree_struct.separator);
  for c=1:Q
    [clpot{c,t}, ll(c)] = normalize_pot(clpot{c,t});
  end
  %logscale(t) = ll(C);
  
  if t >= 2
    phiC = marginalize_pot(clpot{C,t}, engine.interface, engine.maximize);
    phiC = set_domain_pot(phiC, engine.interface+ss); % shift forward to slice 2
    phiD = marginalize_pot(clpot{D,t-1}, engine.interface+ss, engine.maximize);
    ratio = divide_by_pot(phiC, phiD);
    clpot{D,t-1} = multiply_by_pot(clpot{D,t-1}, ratio);
  end
end

loglik = sum(logscale);


%%%%%%%
function [clpot, seppot] = init_pot(cliques, clqs, pots, pot_type, onodes, ns, cnodes);

% Set the clique potentials to all 1s
C = length(cliques);
clpot = cell(1,C);
for i=1:C
  clpot{i} = mk_initial_pot(pot_type, cliques{i}, ns, cnodes, onodes);
end

% Multiply on specified potentials
for i=1:length(clqs)
  c = clqs(i);
  clpot{c} = multiply_by_pot(clpot{c}, pots{i});
end

seppot = cell(C,C); % implicitely initialized to 1


%%%%
function [clpot, seppot] = collect_evidence(clpot, seppot, maximize, postorder, postorder_parents,...
					    separator)
for n=postorder %postorder(1:end-1)
  for p=postorder_parents{n}
    %clpot{p} = divide_by_pot(clpot{n}, seppot{p,n}); % dividing by 1 is redundant
    seppot{p,n} = marginalize_pot(clpot{n}, separator{p,n}, maximize);
    clpot{p} = multiply_by_pot(clpot{p}, seppot{p,n});
  end
end


%%%%
function [clpot, seppot] = distribute_evidence(clpot, seppot, maximize, preorder, preorder_children,...
					       separator)
for n=preorder
  for c=preorder_children{n}
    clpot{c} = divide_by_pot(clpot{c}, seppot{n,c}); 
    seppot{n,c} = marginalize_pot(clpot{n}, separator{n,c}, maximize);
    clpot{c} = multiply_by_pot(clpot{c}, seppot{n,c});
  end
end
