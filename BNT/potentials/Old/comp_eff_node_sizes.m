function ens = comp_eff_node_sizes(ns, cnodes, ev, domain)

dnodes = mysetdiff(1:length(ns), cnodes);
odom = domain(~isemptycell(evidence(domain)));
cdom = myintersect(cnodes, domain);
ddom = myintersect(dnodes, domain);
cobs = myintersect(cdom, odom);
dobs = myintersect(ddom, odom);
ens = ns; 
ens(cobs) = 0;
ens(dobs) = 1;
