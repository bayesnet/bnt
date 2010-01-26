function L = log_marg_prob_node_case(CPD, y, x)
% LOG_MARG_PROB_NODE_CASE Compute prod_m log P(x(i,m)| x(pi_i,m)) for node i (tabular)
% L = log_marg_prob_node_case(CPD, self_ev, parent_ev)
% 
% This is a slightly optimised version of log_marg_prob_node.
% We assume we have exactly 1 case, i.e., y is a scalar and x is a vector (not a cell array).

sz = CPD.sizes;
nparents = length(sz)-1;

% We assume the CPTs are already set to the mean of the posterior (due to update_params)

switch nparents
 case 0, p = CPD.CPT(y);
 case 1, p = CPD.CPT(x(1), y);
 case 2, p = CPD.CPT(x(1), x(2), y);
 case 3, p = CPD.CPT(x(1), x(2), x(3), y);
 otherwise,
  ind = subv2ind(sz, [x y]);
  p = CPD.CPT(ind);
end
L = log(p);
