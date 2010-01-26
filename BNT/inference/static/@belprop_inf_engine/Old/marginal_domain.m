function marginal = marginal_domain(engine, i)
% MARGINAL_DOMAIN Return the marginal on the specified domain (belprop)
% marginal = marginal_domain(engine, i)

marginal = pot_to_marginal(engine.marginal_domains{i});
