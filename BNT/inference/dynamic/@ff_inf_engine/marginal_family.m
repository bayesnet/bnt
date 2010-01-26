function marginal = marginal_family(engine, i, t)
% MARGINAL_FAMILY Compute the marginal on the specified family (ff)
% marginal = marginal_family(engine, i, t)


if engine.filter
  error('can''t currently use marginal_family when filtering with ff');
end

if nargin < 3, t = 1; end

% The method is similar to the following HMM equation:
% xi(i,j,t) = normalise( alpha(i,t) * transmat(i,j) * obsmat(j,t+1) * beta(j,t+1) )
% where xi(i,j,t) = Pr(Q(t)=i, Q(t+1)=j | y(1:T))

bnet = bnet_from_engine(engine);
ss = length(bnet.intra);

if myismember(i, engine.onodes)
  ps = parents(bnet.dag, i);
  p = ps(1);
  marginal = pot_to_marginal(engine.marginals{ps(1),t});
  fam = ([ps i]) + (t-1)*ss;
elseif t==1
  marginal = pot_to_marginal(engine.marginals{i,t});
  fam = i + (t-1)*ss;
else
  pot = engine.CPDpot{i,t};
  c = engine.obschild(i);
  if c>0
    pot = multiply_by_pot(pot, engine.CPDpot{c,t});
  end
  pot = multiply_by_pot(pot, engine.back{i,t});
  ps = parents(bnet.dag, i+ss);
  for p=ps(:)'
    pot = multiply_by_pot(pot, engine.fwd{p,t-1});
  end
  marginal = pot_to_marginal(normalize_pot(pot));
  fam = ([ps i+ss]) + (t-2)*ss;
end

% we convert the domain to the unrolled numbering system
% so that update_ess extracts the right evidence.
marginal.domain = fam;
