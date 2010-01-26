function pot = set_domain_pot(pot, domain)
% SET_DOMAIN_POT Change the domain of a potential (cgpot)
% pot = set_domain_pot(pot, domain)

delta = domain(1) - pot.domain(1);
assert(all(domain == pot.domain + delta));
pot.domain = pot.domain + delta;
pot.ddom = pot.ddom + delta;
pot.cdom = pot.cdom + delta;
