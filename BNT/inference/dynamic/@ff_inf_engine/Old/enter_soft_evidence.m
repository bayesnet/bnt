function [marginals, loglik] = enter_soft_evidence(engine, CPDpot, observed, pot_type, filter)
% ENTER_SOFT_EVIDENCE Add the specified soft evidence to the network (bk_ff)
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
local_logscale = zeros(1,length(hnodes));

t = 1;
for i=hnodes
  fwd{i,t} = CPDpot{i,t};
end
for i=onodes
  p = parents(bnet.dag, i);
  assert(length(p)==1);
  ev = marginalize_pot(CPDpot{i,t}, p);
  fwd{p,t} = multiply_by_pot(fwd{p,t}, ev);
end
for i=hnodes
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
  end
  for i=onodes
    p = parents(bnet.dag, i);
    assert(length(p)==1);
    temp = pot_to_marginal(CPDpot{i,t}); 
    ev = dpot(p, ns(p), temp.T);
    fwd{p,t} = multiply_by_pot(fwd{p,t}, ev);
  end
  
  for i=hnodes
    [fwd{i,t}, local_logscale(i)] = normalize_pot(fwd{i,t});
  end
  logscale(t) = sum(local_logscale);
end

marginals = fwd;
loglik = sum(logscale);
