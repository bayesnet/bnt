function pot = set_domain_pot(pot, domain)
% SET_DOMAIN_POT Change the domain of a potential (cgpot)
% pot = set_domain_pot(pot, domain)

delta = domain(1) - pot.domain(1);
assert(all(domain == pot.domain + delta));
pot.domain = pot.domain + delta;
pot.ddom = pot.ddom + delta;
pot.cdom = pot.cdom + delta;
cdomain = pot.cdom;
n = prod(pot.dsizes);
if(pot.subtype == 'm')
    for i = 1: n
        pot.mom{i} = set_domain_pot(pot.mom{i}, cdomain);
    end
end
if(pot.subtype == 'c')
    for i = 1: n
        pot.can{i} = set_domain_pot(pot.can{i}, cdomain);
    end
end
