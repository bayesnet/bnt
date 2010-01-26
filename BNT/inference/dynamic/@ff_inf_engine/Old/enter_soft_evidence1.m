function [marginals, loglik] = enter_soft_evidence(engine, CPDpot, observed, pot_type, filter)
% ENTER_SOFT_EVIDENCE Add the specified soft evidence to the network (ff)
% [marginals, loglik] = enter_soft_evidence(engine, CPDpot, observed, pot_type, filter)

assert(pot_type == 'd');
[ss T] = size(CPDpot);
fwd = cell(ss,T);
hnodes = engine.hnodes(:)';
onodes = engine.onodes(:)';
bnet = bnet_from_engine(engine);
ns = bnet.node_sizes;
onodes2 = [onodes onodes+ss];
ns(onodes2) = 1;

logscale = zeros(1,T);
H = length(hnodes);
local_logscale = zeros(1,ss);

obschild = zeros(1,ss);
for i=hnodes
  ocs = myintersect(children(bnet.dag, i), onodes);
  assert(length(ocs)==1);
  obschild(i) = ocs(1);
end  
  
t = 1;
for i=hnodes
  fwd{i,t} = CPDpot{i,t};
  c = obschild(i);
  temp = pot_to_marginal(CPDpot{c,t}); 
  ev = dpot(i, ns(i), temp.T);
  fwd{i,t} = multiply_by_pot(fwd{i,t}, ev);
  [fwd{i,t}, local_logscale(i)] = normalize_pot(fwd{i,t});
end
logscale(t) = sum(local_logscale);

for t=2:T
  for i=hnodes
    ps = parents(bnet.dag, i+ss);
    assert(all(ps<=ss)); % in previous slice
    prior = CPDpot{i,t};
    for p=ps(:)'
      prior = multiply_by_pot(prior, fwd{p,t-1});
    end
    fwd{i,t} = marginalize_pot(prior, i+ss);
    fwd{i,t} = set_domain_pot(fwd{i,t}, i);
    c = obschild(i);
    temp = pot_to_marginal(CPDpot{c,t});
    ev = dpot(i, ns(i), temp.T);
    fwd{i,t} = multiply_by_pot(fwd{i,t}, ev);
    [fwd{i,t}, local_logscale(i)] = normalize_pot(fwd{i,t});
  end
  logscale(t) = sum(local_logscale);
end

loglik = sum(logscale);


if filter
  marginals = fwd;
  return;
end

back = cell(ss,T);
t = T;
for i=hnodes
  back{i,t} = dpot(i, ns(i));
  back{i,t} = set_domain_pot(back{i,t}, i+ss);
end
for t=T-1:-1:1
  for i=hnodes
    pot = CPDpot{i,t+1};
    pot = multiply_by_pot(pot, back{i,t+1});
    c = obschild(i);
    temp = pot_to_marginal(CPDpot{c,t+1});
    ev = dpot(i, ns(i), temp.T);
    pot = multiply_by_pot(pot, ev);
    back{i,t} = marginalize_pot(pot, i);
    back{i,t} = normalize_pot(back{i,t});
    back{i,t} = set_domain_pot(back{i,t}, i+ss);
  end
end



% COMBINE
for t=1:T
  for i=hnodes
    back{i,t} = set_domain_pot(back{i,t}, i);
    fwd{i,t} = multiply_by_pot(fwd{i,t}, back{i,t});
    marginals{i,t} = normalize_pot(fwd{i,t});
    %fwdback{i,t} = normalize_pot(multiply_pots(fwd{i,t}, back{i,t}));
  end
end
