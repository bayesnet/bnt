function pot = mk_initial_pot(pot_type, dom, ns, cnodes, onodes)
% MK_INITIAL_POT A "initial" potential is one which has not had any evidence entered into it.
% pot = mk_initial_pot(pot_type, domain, node_sizes, cnodes, onodes)
%
% pot_type is one of 'd', 'g', 'cg' or 'u'
% domain is the set of nodes to be included in the potential.
% node_sizes(i) is the size of node i.

switch pot_type 
 case 'd',
  ns(onodes) = 1;
  pot = dpot(dom, ns(dom));
 case 'u',
  ns(onodes) = 1;
  pot = upot(dom, ns(dom));
 case 'g',
  ns(onodes) = 0;
  pot = cpot(dom, ns(dom));
 case 'cg',
  dnodes = mysetdiff(1:length(ns), cnodes);
  ddom = myintersect(dnodes, dom);
  cdom = myintersect(cnodes, dom);
  dobs = myintersect(dnodes, onodes);
  cobs = myintersect(cnodes, onodes);
  ns(dobs) = 1;
  ns(cobs) = 0;
  pot = cgpot(ddom, cdom, ns);
end

