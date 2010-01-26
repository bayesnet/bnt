function marginal = marginal_family(engine, i, t)
% MARGINAL_FAMILY Compute the marginal on the specified family (jtree_unrolled_dbn)
% marginal = marginal_family(engine, i, t)

if nargin < 3, t = 1; end

bnet = bnet_from_engine(engine);
ss = length(bnet.intra);
marginal = marginal_family(engine.sub_engine, i + (t-1)*ss);
 
