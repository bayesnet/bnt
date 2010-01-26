function pi = CPD_to_pi(CPD, msg_type, n, ps, msg, evidence)
% CPD_TO_PI Compute pi vector (noisyor)
% pi = CPD_to_pi(CPD, msg_type, n, ps, msg)
% Pearl p188 eqn 4.57
  
switch msg_type
 case 'd',
   pi = sum_prod_CPD_and_pi_msgs(CPD, n, ps, msg);
   pi = [pi 1-pi]';
 case 'g', 
  error('can''t convert noisy-or CPD to Gaussian pi')
end
