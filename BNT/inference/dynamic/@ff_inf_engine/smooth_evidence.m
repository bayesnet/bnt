function [marginals, fwd, back, loglik] = smooth_evidence(engine, CPDpot, observed, pot_type)
% [marginals, fwd, back, loglik] = smooth_evidence(engine, CPDpot, observed, pot_type) (ff)

error('ff smoothing is broken');

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
  
t = 1;
for i=hnodes
  fwd{i,t} = CPDpot{i,t};
  c = engine.obschild(i);
  if 0 %  c > 0
    fwd{i,t} = multiply_by_pot(fwd{i,t}, CPDpot{c, t});
  end
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
    c = engine.obschild(i);
    if 0 % c > 0
      fwd{i,t} = multiply_by_pot(fwd{i,t}, CPDpot{c,t});
    end
    [fwd{i,t}, local_logscale(i)] = normalize_pot(fwd{i,t});
  end
  logscale(t) = sum(local_logscale);
end

loglik = sum(logscale);

back = cell(ss,T);
t = T;
for i=hnodes
  pot = dpot(i, ns(i));
  cs = children(bnet.intra, i);
  for c=cs(:)'
    pot = multiply_pots(pot, CPDpot{c,t});
  end
  back{i,t} = marginalize_pot(pot, i);
  back{i,t} = normalize_pot(back{i,t});
  back{i,t} = set_domain_pot(back{i,t}, i+ss);
end
for t=T-1:-1:1
  for i=hnodes
    pot = dpot(i, ns(i));
    cs = children(bnet.inter, i);
    for c=cs(:)'
      pot = multiply_pots(pot, back{c,t+1});
      pot = multiply_pots(pot, CPDpot{c,t+1});
    end
    cs = children(bnet.intra, i);
    for c=cs(:)'
      pot = multiply_pots(pot, CPDpot{c,t});
    end
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
  end
end
