function [policy, EU] = upot_to_opt_policy(pot)
% UPOT_TO_OPT_POLICY Compute an optimal deterministic policy given a utility potential
% [policy, EU] = upot_to_opt_policy(pot)
%
% policy(a,b, ..., z) = P(do z | a, b, ..), which will be a delta function
% EU is the contraction of this potential, i.e., P .* U

sz = pot.sizes; % mysize(pot.p);
if isempty(sz)
  EU = pot.u;
  policy = [];
  return;
end

parent_size = prod(sz(1:end-1));
self_size = sz(end); 
C = pot.p .* pot.u; % contraction
C = reshape(C, parent_size, self_size);
policy = zeros(parent_size, self_size);
for i=1:parent_size
  act = argmax(C(i,:));
  policy(i, act) = 1;
end
policy = myreshape(policy, sz);
EU = sum(C(:));
