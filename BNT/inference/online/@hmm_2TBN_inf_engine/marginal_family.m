function marginal = marginal_family(engine, b, i, t, add_ev)
% MARGINAL_FAMILY Compute the marginal on the specified family (hmm_2TBN)
% marginal = marginal_family(engine, b, i, t, add_ev)

ns = engine.eff_node_sizes(:);
ss = engine.slice_size;

if t==1 % | ~engine.persist_bitv(i)
  bigT = b.gamma;
  ps = engine.parents{i};
  dom = [ps i];
  %id = engine.marg_fam1_ndx_id(i);
  bigdom = 1:ss;
  bigsz = ns(bigdom);
  bigdom = bigdom + (t-1)*ss;
else % some parents are in previous slice
  bigT = b.xi; % (t-1,t)
  ps = engine.parents{i+ss};
  dom = [ps i+ss] + (t-2)*ss;
  %id = engine.marg_fam2_ndx_id(i);
  bigdom = 1:(2*ss); % domain of xi(:,:,t)
  bigsz = ns(bigdom);
  bigdom = bigdom + (t-2)*ss;
end
marginal.domain = dom;

%ndx = get_ndx(id, engine.ndx_type);
%marginal.T = marg_table_ndx(bigT, engine.maximize, ndx, engine.ndx_type);
%global SD_NDX
%ndx = SD_NDX{id};
%marginal.T = marg_table_ndxSD(bigT, engine.maximize, ndx);
marginal.T = marg_table(bigT, bigdom, bigsz, dom, engine.maximize); 

assert(~add_ev);

