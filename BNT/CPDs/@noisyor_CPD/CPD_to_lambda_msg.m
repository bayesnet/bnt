function lam_msg = CPD_to_lambda_msg(CPD, msg_type, n, ps, msg, p, evidence)
% CPD_TO_LAMBDA_MSG Compute lambda message (noisyor)
% lam_msg = CPD_to_lambda_msg(CPD, msg_type, n, ps, msg, p)
% Pearl p190 top eqn

switch msg_type
  case 'd', 
   l0 = msg{n}.lambda(1);
   l1 = msg{n}.lambda(2);
   Pi = sum_prod_CPD_and_pi_msgs(CPD, n, ps, msg, p);
   i = find(p==ps); % p is n's i'th parent
   q = CPD.inhibit(i);
   lam_msg = zeros(2,1);
   for u=0:1
     lam_msg(u+1) = l1 - (q^u)*(l1 - l0)*Pi;
   end       
 case 'g',
  error('noisyor_CPD can''t create Gaussian msgs')
end
