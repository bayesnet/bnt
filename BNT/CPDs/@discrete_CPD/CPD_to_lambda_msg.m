function lam_msg = CPD_to_lambda_msg(CPD, msg_type, n, ps, msg, p, evidence)
% CPD_TO_LAMBDA_MSG Compute lambda message (discrete)
% lam_msg = compute_lambda_msg(CPD, msg_type, n, ps, msg, p, evidence)
% Pearl p183 eq 4.52

switch msg_type
  case 'd',
   T = prod_CPT_and_pi_msgs(CPD, n, ps, msg, p);
   mysize = length(msg{n}.lambda);
   lambda = dpot(n, mysize, msg{n}.lambda);
   T = multiply_by_pot(T, lambda);
   lam_msg = pot_to_marginal(marginalize_pot(T, p));
   lam_msg = lam_msg.T;           
 case 'g',
  error('discrete_CPD can''t create Gaussian msgs')
end
