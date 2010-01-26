function m = marginal_family(engine, n, add_ev)
% MARGINAL_FAMILY Compute the marginal on i's family (loopy)
% m = marginal_family(engine, n, add_ev)

if nargin < 3, add_ev = 0; end

bnet = bnet_from_engine(engine);
ns = bnet.node_sizes;
ps = parents(bnet.dag, n);
dom = [ps n];
CPD = bnet.CPD{bnet.equiv_class(n)};

switch engine.msg_type
  case 'd',
   % The method is similar to the following HMM equation:
   % xi(i,j,t) = normalise( alpha(i,t) * transmat(i,j) * obsmat(j,t+1) * beta(j,t+1) )
   % where xi(i,j,t) = Pr(Q(t)=i, Q(t+1)=j | y(1:T))   
   % beta == lambda, alpha == pi, alpha from each parent = pi msg
   % In general, if A,B are parents of C,
   % P(A,B,C) = P(C|A,B) pi_msg(A->C) pi_msg(B->C) lambda(C)
   % where lambda(C) = P(ev below and including C|C) = prod incoming lamba_msg(children->C)
   % and pi_msg(X->C) = P(X|ev above) etc
   
   T = dpot(dom, ns(dom), CPD_to_CPT(CPD));
   for j=1:length(ps)
     p = ps(j);
     pi_msg = dpot(p, ns(p), engine.msg{n}.pi_from_parent{j});
     T = multiply_by_pot(T, pi_msg);
   end         
   lambda = dpot(n, ns(n), engine.msg{n}.lambda);
   T = multiply_by_pot(T, lambda);
   T = normalize_pot(T);
   m = pot_to_marginal(T);
   if ~add_ev
     m.T = shrink_obs_dims_in_table(m.T, dom, engine.evidence);
   end
 case 'g',
  if engine.disconnected_nodes_bitv(n)
    m.T = 1;
    m.domain = dom;
    if add_ev
      m = add_ev_to_dmarginal(m, engine.evidence, ns)
    end
    return;
  end

  [m, C, W] = gaussian_CPD_params_given_dps(CPD, dom, engine.evidence);
  cdom = myintersect(dom, bnet.cnodes);
  pot = linear_gaussian_to_cpot(m, C, W, dom, ns, cdom, engine.evidence); 
  % linear_gaussian_to_cpot will set the effective size of observed nodes to 0,
  % so we need to do this explicitely for the messages, too,
  % so they are all the same size.
  obs_bitv = ~isemptycell(engine.evidence);
  ps = parents(engine.msg_dag, n);
  for j=1:length(ps)
    p = ps(j);
    msg = engine.msg{n}.pi_from_parent{j};
    if obs_bitv(p)
      pi_msg = mpot(p, 0);
    else
      pi_msg = mpot(p, ns(p), 0, msg.mu, msg.Sigma);
    end
    pot = multiply_by_pot(pot, mpot_to_cpot(pi_msg));
  end         
  msg = engine.msg{n}.lambda;
  if obs_bitv(n)
    lambda = cpot(n, 0);
  else
    lambda = cpot(n, ns(n), 0, msg.info_state, msg.precision);
  end
  pot = multiply_by_pot(pot, lambda);
  m = pot_to_marginal(pot);
  if add_ev
    m = add_evidence_to_gmarginal(m, engine.evidence, bnet.node_sizes, bnet.cnodes);
  end
end




