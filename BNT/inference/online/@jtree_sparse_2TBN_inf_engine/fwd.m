function [f, logscale] = fwd(engine, fpast, ev, t)
% Forwards pass.

bnet = bnet_from_engine(engine);
ss = bnet.nnodes_per_slice;

ev2 = cell(ss, 2);
ev2(:,1) = fpast.evidence;
ev2(:,2) = ev;
CPDpot = cell(1,ss);
for n=1:ss
  fam = family(bnet.dag, n, 2);
  e = bnet.equiv_class(n, 2);
  CPDpot{n} = convert_to_pot(bnet.CPD{e}, engine.pot_type, fam(:), ev2);
end       
f.evidence = ev;
f.t = t;

% get prior
int = engine.interface;
if fpast.t==1
  prior = marginalize_pot(fpast.clpot{engine.int_clq1}, int, engine.maximize);
else
  prior = marginalize_pot(fpast.clpot{engine.out_clq}, int+ss, engine.maximize);
  prior = set_domain_pot(prior, int); % shift back to slice 1
end

pots = [ {prior} CPDpot ];
slice1 = 1:ss;
slice2 = slice1 + ss; 
CPDclqs = engine.clq_ass_to_node(slice2);
D = engine.in_clq;
clqs = [D CPDclqs];

[f.clpot, f.seppot] =  init_pot(engine.jtree_engine, clqs, pots, engine.pot_type, engine.observed);
[f.clpot, f.seppot] = collect_evidence(engine.jtree_engine, f.clpot, f.seppot);
for c=1:length(f.clpot)
  if isa(f.clpot{c}, 'struct')
     domain = f.clpot{c}.domain;
     sizes = f.clpot{c}.sizes;
     T = f.clpot{c}.T;
     f.clpot{c} = dpot(domain, sizes, T);
  end
  [f.clpot{c}, ll(c)] = normalize_pot(f.clpot{c});
end
logscale = ll(engine.root);

