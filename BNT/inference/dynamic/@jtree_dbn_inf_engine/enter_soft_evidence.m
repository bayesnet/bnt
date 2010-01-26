function [clpot, loglik] = enter_soft_evidence(engine, CPDpot, observed, pot_type)
% ENTER_SOFT_EVIDENCE Add the specified soft evidence to the network (jtree_dbn)
% [clpot, loglik] = enter_soft_evidence(engine, CPDpot, observed, pot_type)

scale = 1;
verbose = 0;

[ss T] = size(CPDpot);
Q = length(engine.jtree_struct.cliques);
clpot = cell(Q,T); % clpot{t} contains evidence from slices (t-1, t) 
seppot = cell(Q,Q,T);
ll = zeros(1,Q);
logscale = zeros(1,T);
bnet = bnet_from_engine(engine);
root = engine.jtree_struct.root_clq;

% Forwards pass.
% Compute distribution on clq C,
% where C is the out interface to (t-1,t).
% Then pass this to clq D, where D is the in inferface to (t+1,t).

% Then propagate from D to later slices.

slice1 = 1:ss;
slice2 = slice1 + ss; 
transient = engine.transient;
persist = engine.persist;
Ntransient = length(transient);
trans = cell(Ntransient,1);
if verbose, fprintf('forward pass\n'); end
for t=1:T
  if verbose, fprintf('%d ', t); end
  if t==1
    pots = [CPDpot(:,1); CPDpot(persist, 2)];
    clqs = engine.jtree_struct.clq_ass_to_node([slice1 persist+ss]);
    obs = find(observed(:,1:2));
  elseif t==T
    clqs = [engine.in_clq1 engine.jtree_struct1.clq_ass_to_node(transient)];
    phi = set_domain_pot(phi, engine.interface); % shift back to slice 1
    for i=1:Ntransient
      trans{i} = CPDpot{transient(i), t};
      trans{i} = set_domain_pot(trans{i}, domain_pot(trans{i})-ss); % shift back to slice 1
    end
    pots = [ {phi}; trans]; 
    obs = find(observed(:,T));
  else
    clqs = [engine.in_clq engine.jtree_struct.clq_ass_to_node([transient persist+ss])];
    phi = set_domain_pot(phi, engine.interface); % shift back to slice 1
    for i=1:Ntransient
      trans{i} = CPDpot{transient(i), t};
      trans{i} = set_domain_pot(trans{i}, domain_pot(trans{i})-ss); % shift back to slice 1
    end
    pots = [ {phi}; trans; CPDpot(persist, t+1)]; 
    obs = find(observed(:,t:t+1));
  end

  if t < T
    [clpot(1:Q,t), seppot(1:Q,1:Q,t)] =  init_pot(engine.jtree_engine, clqs, pots, pot_type, obs);
    [clpot(1:Q,t), seppot(1:Q,1:Q,t)] = collect_evidence(engine.jtree_engine, clpot(1:Q,t), seppot(1:Q,1:Q,t));
  else
    Q = length(engine.jtree_struct1.cliques);
    root = engine.jtree_struct1.root_clq;
    [clpot(1:Q,t), seppot(1:Q,1:Q,t)] =  init_pot(engine.jtree_engine1, clqs, pots, pot_type, obs);
    [clpot(1:Q,t), seppot(1:Q,1:Q,t)] = collect_evidence(engine.jtree_engine1, clpot(1:Q,t), seppot(1:Q,1:Q,t));
  end


  if scale
  for c=1:Q
    [clpot{c,t}, ll(c)] = normalize_pot(clpot{c,t});
  end
  logscale(t) = ll(root);
  end
  
  if t < T
    % bug fix by Bob Welch 30 Jan 04
    phi = marginalize_pot(clpot{engine.out_clq,t}, engine.interface+ss,engine.maximize);
    %phi = marginalize_pot(clpot{root,t}, engine.interface+ss, engine.maximize);
  end
end

if scale
loglik = sum(logscale);
else
loglik = [];
end


% Backwards pass.
% Pass evidence from clq C to clq D,
% where C is the in interface to (t,t+1) and D is the out inferface to (t-1,t)
% Then propagate evidence from D to earlier slices.
% (C and D are reversed names from the tech report!)
D = engine.out_clq;
if verbose, fprintf('\nbackwards pass\n'); end
for t=T:-1:1
  if verbose, fprintf('%d ', t); end
  
  if t == T
    Q = length(engine.jtree_struct1.cliques);
    C = engine.in_clq1;
    [clpot(1:Q,t), seppot(1:Q,1:Q,t)] = distribute_evidence(engine.jtree_engine1, clpot(1:Q,t), seppot(1:Q,1:Q,t));
  else
    Q = length(engine.jtree_struct.cliques);
    C = engine.in_clq;
    [clpot(1:Q,t), seppot(1:Q,1:Q,t)] = distribute_evidence(engine.jtree_engine, clpot(1:Q,t), seppot(1:Q,1:Q,t));
  end

  if scale
  for c=1:Q
    [clpot{c,t}, ll(c)] = normalize_pot(clpot{c,t});
  end
  end
  
  if t >= 2
    phiC = marginalize_pot(clpot{C,t}, engine.interface, engine.maximize);
    phiC = set_domain_pot(phiC, engine.interface+ss); % shift forward to slice 2
    phiD = marginalize_pot(clpot{D,t-1}, engine.interface+ss, engine.maximize);
    ratio = divide_by_pot(phiC, phiD);
    clpot{D,t-1} = multiply_by_pot(clpot{D,t-1}, ratio);
  end
end
if verbose, fprintf('\n'); end



