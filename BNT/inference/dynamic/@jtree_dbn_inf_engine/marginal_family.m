function m = marginal_family(engine, i, t, add_ev)
% MARGINAL_FAMILY Compute the marginal on the specified family (jtree_dbn)
% marginal = marginal_family(engine, i, t)

% This is just like inf_engine/marginal_family, except when we call
% marginal_nodes, we provide a 4th argument, to tell it's a family.

if nargin < 3, t = 1; end
if nargin < 4, add_ev = 0; end

bnet = bnet_from_engine(engine);
if t==1
  m = marginal_nodes(engine, family(bnet.dag, i), t, add_ev, 1);
else
  ss = length(bnet.intra);
  fam = family(bnet.dag, i+ss);
  if any(fam<=ss)
    % i has a parent in the preceeding slice
    % Hence the lowest numbered slice containing the family is t-1
    m = marginal_nodes(engine, fam, t-1, add_ev, 1);
  else
    % The family all fits inside slice t
    % Hence shift the indexes back to slice 1
    m = marginal_nodes(engine, fam-ss, t, add_ev, 1);
  end
end     
