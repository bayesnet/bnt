function marginal = marginal_family(engine, i, t, add_ev)
% MARGINAL_FAMILY Compute the marginal on the specified family (hmm)
% marginal = marginal_family(engine, i, t, add_ev)

if nargin < 3, t = 1; end
if nargin < 4, add_ev = 0; end

ns = engine.eff_node_sizes(:);
ss = engine.slice_size;

if t==1 | ~engine.persist_bitv(i)
  bigT = engine.one_slice_marginal(:,t);
  ps = engine.parents{i};
  dom = [ps i] + (t-1)*ss;
  bigdom = 1:ss;
  bigsz = ns(bigdom);
  bigdom = bigdom + (t-1)*ss;
else % some parents are in previous slice
  bigT = engine.two_slice_marginal(:,t-1); % t-1 and t
  ps = engine.parents{i+ss};
  dom = [ps i+ss] + (t-2)*ss; 
  bigdom = 1:(2*ss); % domain of xi(:,:,t)
  bigsz = ns(bigdom);
  bigdom = bigdom + (t-2)*ss;
end
marginal.domain = dom;

marginal.T = marg_table(bigT, bigdom, bigsz, dom, engine.maximize); 
marginal.mu = []; 
marginal.Sigma = [];

if add_ev
  marginal = add_ev_to_dmarginal(marginal, engine.evidence, engine.node_sizes);
end    

