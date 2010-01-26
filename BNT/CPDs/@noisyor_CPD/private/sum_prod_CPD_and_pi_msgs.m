function pi = sum_prod_CPD_and_pi_msgs(CPD, n, ps, msg, except)
% SUM_PROD_CPD_AND_PI_MSGS Compute pi = sum_{u\p} P(n|u) prod_{ui in ps\p} pi_msg(ui->n)
% pi = sum_prod_CPD_and_pi_msgs(CPD, n, ps, msg, p)
%
% pi  = prod_i (qi pi_msg(ui->n) + 1 - pi_msg(ui->n)) = prod_i (1 - ci pi_msg(ui->n))
% is the product of the endorsement withheld (Pearl p188 eqn 4.56)
% We skip p from this product, if specified.

if nargin < 5, except = -1; end
pi = 1;
for i=1:length(ps)
  p = ps(i);
  if p ~= except
    pi_from_parent = msg{n}.pi_from_parent{i};
    q = CPD.inhibit(i);
    c = 1-q;
    pi = pi * (1 - c*pi_from_parent(2));
  end
end
% The pi msg that a leak node sends to its child is [0 1]
% since its own pi is [0 1] and its lambda to self is [0 1].
q = CPD.leak_inhibit;
% 1 - c*pi_from_parent = 1-c*1 = q
pi = pi * q;
                 
