function pi = CPD_to_pi(CPD, msg_type, n, ps, msg, evidence)
% COMPUTE_PI Compute pi vector (discrete) 
% pi = compute_pi(CPD, msg_type, n, ps, msg, evidence)
% Pearl p183 eq 4.51

switch msg_type
  case 'd',
   T = prod_CPT_and_pi_msgs(CPD, n, ps, msg);
   pi = pot_to_marginal(marginalize_pot(T, n));
   pi = pi.T(:);                   
 case 'g', 
  error('can only convert discrete CPD to Gaussian pi if observed')
end
