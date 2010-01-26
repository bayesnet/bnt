function [fwdback, loglik, fwd_frontier, back_frontier] = enter_soft_evidence(engine, CPD, onodes, pot_type, filter)
% ENTER_SOFT_EVIDENCE Add soft evidence to network (frontier)
% [fwdback, loglik] = enter_soft_evidence(engine, CPDpot, onodes, filter)

if nargin < 3, filter = 0; end

[ss T] = size(CPD);
bnet = bnet_from_engine(engine);
ns = repmat(bnet.node_sizes_slice(:), 1, T);
cnodes = unroll_set(bnet.cnodes(:), ss, T);

% FORWARDS
fwd = cell(ss,T);
ll = zeros(1,T);
S = 2*ss; % num. intermediate frontiers to get from t to t+1
frontier = cell(S,T);

% Start with empty frontier, and add each node in slice 1
init = mk_initial_pot(pot_type, [], ns, cnodes, onodes);  
t = 1;
s = 1;
j = 1;
frontier{s,t} = update(init, j, 1, CPD{j}, engine.fdom1{s}, pot_type, ns, cnodes, onodes);
fwd{j} = frontier{s,t};
for s=2:ss
  j = s; % add node j at step s
  frontier{s,t} = update(frontier{s-1,t}, j, 1, CPD{j}, engine.fdom1{s}, pot_type, ns, cnodes, onodes);
  fwd{j} = frontier{s,t};
end
frontier{S,t} = frontier{ss,t};
[frontier{S,t}, ll(1)] = normalize_pot(frontier{S,t});

% Now move frontier from slice to slice
OPS = engine.ops;
add = OPS>0;
nodes = [zeros(S,1) unroll_set(abs(OPS(:)), ss, T-1)];
for t=2:T
  offset = (t-2)*ss;
  for s=1:S
    if s==1
      prev_ndx = (t-2)*S + S; % S,t-1
    else
      prev_ndx = (t-1)*S + s-1; % s-1,t
    end
    j = nodes(s,t);
    frontier{s,t} = update(frontier{prev_ndx}, j, add(s), CPD{j}, engine.fdom{s}+offset, pot_type, ns, cnodes, onodes);
    if add(s)
      fwd{j} = frontier{s,t};
    end
  end
  [frontier{S,t}, ll(t)] = normalize_pot(frontier{S,t});
end
loglik = sum(ll);


fwd_frontier = frontier;

if filter
  fwdback = fwd;
  return;
end


% BACKWARDS
back = cell(ss,T);
add = ~add; % forwards add = backwards remove 
frontier = cell(S,T+1);
t = T;
dom = (1:ss) + (t-1)*ss;
frontier{1,T+1} = mk_initial_pot(pot_type, dom, ns, cnodes, onodes); % all 1s for last slice
for t=T:-1:2
  offset = (t-2)*ss;
  for s=S:-1:1 % reverse order
    if s==S
      prev_ndx = t*S + 1; % 1,t+1
    else
      prev_ndx = (t-1)*S + (s+1); % s+1,t
    end
    j = nodes(s,t);
    if ~add(s)
      back{j} = frontier{prev_ndx}; % save frontier before removing
    end
    frontier{s,t} = rev_update(frontier{prev_ndx}, t, s, j, add(s), CPD{j}, engine.fdom{s}+offset, pot_type, ns, cnodes, onodes);
  end
  frontier{1,t} = normalize_pot(frontier{1,t});
end
% Remove each node in first slice until left with empty set
t = 1;
frontier{ss+1,t} = frontier{1,2};
add = 0;
for s=ss:-1:1
  j = s; % remove node j at step s
  back{j} = frontier{s+1,t};
  frontier{s,t} = rev_update(frontier{s+1,t}, t, s, j, add, CPD{j}, 1:s, pot_type, ns, cnodes, onodes);
end

% COMBINE
for t=1:T
  for i=1:ss
    %fwd{i,t} = multiply_by_pot(fwd{i,t}, back{i,t});
    %fwdback{i,t} = normalize_pot(fwd{i,t});
    fwdback{i,t} = normalize_pot(multiply_pots(fwd{i,t}, back{i,t}));
  end
end

back_frontier = frontier;

%%%%%%%%%%
function new_frontier = update(old_frontier, j, add, CPD, newdom, pot_type, ns, cnodes, onodes)

if add
  new_frontier = mk_initial_pot(pot_type, newdom, ns, cnodes, onodes);      
  new_frontier = multiply_by_pot(new_frontier, old_frontier);
  new_frontier = multiply_by_pot(new_frontier, CPD);
else
  new_frontier = marginalize_pot(old_frontier, mysetdiff(domain_pot(old_frontier), j));    
end


%%%%%%
function new_frontier = rev_update(old_frontier, t, s, j, add, CPD, junk, pot_type, ns, cnodes, onodes)

olddom = domain_pot(old_frontier);
assert(isequal(junk, olddom));

if add
  % add: extend domain to include j by multiplying by 1
  newdom = myunion(olddom, j);
  new_frontier = mk_initial_pot(pot_type, newdom, ns, cnodes, onodes);      
  new_frontier = multiply_by_pot(new_frontier, old_frontier);
  %fprintf('t=%d, s=%d, add %d to %s to make %s\n', t, s, j, num2str(olddom), num2str(newdom));
else 
  % remove: multiply in CPT and then marginalize out j
  % parents of j are guaranteed to be in old_frontier, else couldn't have added j on fwds pass
  old_frontier = multiply_by_pot(old_frontier, CPD);
  newdom = mysetdiff(olddom, j);
  new_frontier = marginalize_pot(old_frontier, newdom);
  %newdom2 = domain_pot(new_frontier);
  %fprintf('t=%d, s=%d, rem %d from %s to make %s\n', t, s, j, num2str(olddom), num2str(newdom2));
end

       
