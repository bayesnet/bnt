function marginal = marginal_family(engine, i, t, add_ev)
% MARGINAL_FAMILY Compute the marginal on the specified family  (hmm)
% marginal = marginal_nodes(engine, i, t, add_ev)
%

if nargin < 3, t = 1; end
if nargin < 4, add_ev = 0; end

bnet = bnet_from_engine(engine);
ss = length(bnet.intra);
if t==1
 fam = family(bnet.dag, i);
 bigpot = engine.one_slice_marginal{t};
 nodes = fam;
else
  fam = family(bnet.dag, i+ss);
  if any(fam <= ss) % family spans 2 slices
    bigpot = engine.two_slice_marginal{t-1}; % t-1 and t
    nodes = fam + (t-2)*ss;
  else
    bigpot = engine.one_slice_marginal{t};
    nodes = fam-ss + (t-1)*ss;
  end
end

marginal = pot_to_marginal(marginalize_pot(bigpot, nodes, engine.maximize));

if add_ev
  marginal = add_ev_to_dmarginal(marginal, engine.evidence, engine.node_sizes);
end    

