function CPT = noisyORtoTable(inhibit, leak_inhibit)
% NOISYORTOTABLE Convert noisyOR distribution to CPT
% function CPT = noisyORtoTable(inhibit, leak_inhibit)
%
% inhibit(i) = prob i'th parent will be inhibited (flipped from 1 to 0)
% leak_inhibit - optional suppression of leak
% CPT(U1,...,Un, X) = Pr(X|U1,...,Un) where the Us are the parents (excluding leak).
% State 1 = off, 2 = on

if nargin < 2, leak_inhibit = 1; end

q = [leak_inhibit inhibit(:)'];

if length(q)==1
  CPT = [q  1-q];
  return;
end

n = length(q);
Bn = ind2subv(2*ones(1,n), 1:(2^n))-1;  % all n bit vectors, with the left most column toggling fastest (LSB)
CPT = zeros(2^n, 2);
% Pr(X=0 | U_1 .. U_n) = prod_{i: U_i = on} q_i =  prod_i q_i ^ U_i = exp(u' * log(q_i))
% This method is problematic when q contains zeros

Q = repmat(q(:)', 2^n, 1);
Q(logical(~Bn)) = 1;
CPT(:,1) = prod(Q,2);
CPT(:,2) = 1-CPT(:,1);

CPT = reshape(CPT(2:2:end), 2*ones(1,n)); % skip cases in which the leak is off       


