function marginal = marginal_family(engine, i, t)
% MARGINAL_FAMILY Compute the marginal on the specified family (ff)
% marginal = marginal_family(engine, i, t)

if nargin < 3, t = 1; end

% The method is similar to the following HMM equation:
% xi(i,j,t) = normalise( alpha(i,t) * transmat(i,j) * obsmat(j,t+1) * beta(j,t+1) )
% where xi(i,j,t) = Pr(Q(t)=i, Q(t+1)=j | y(1:T))

bnet = bnet_from_engine(engine);

if myismember(i, engine.onodes)
  ps = parents(bnet.dag, i);
  p = ps(1);
  marginal = pot_to_marginal(engine.marginals{p,t});
  marginal.domain = [p i];
  return;
end

if t==1
  marginal = pot_to_marginal(engine.marginals{i,t});
  return;
end

bnet = bnet_from_engine(engine);
ss = length(bnet.intra);
pot = engine.CPDpot{i,t};
c = engine.obschild(i);
pot = multiply_by_pot(pot, engine.CPDpot{c,t});
pot = multiply_by_pot(pot, engine.back{i,t});
ps = parents(bnet.dag, i+ss);
for p=ps(:)'
  pot = multiply_by_pot(pot, engine.fwd{p,t-1});
end
marginal = pot_to_marginal(normalize_pot(pot));


