function smallpot = marginalize_pot(bigpot, keep, maximize, useC)
% MARGINALIZE_POT Marginalize a cpot onto a smaller domain.
% smallpot = marginalize_pot(bigpot, keep, maximize, useC)
%
% The maximize argument is ignored - maxing out a Gaussian is the same as summing it out,
% since the mode and mean are equal.
% The useC argument is ignored.

node_sizes = sparse(1, max(bigpot.domain));
node_sizes(bigpot.domain) = bigpot.sizes;
sum_over = mysetdiff(bigpot.domain, keep);

if sum(node_sizes(sum_over))==0 % isempty(sum_over)
  %smallpot = bigpot;
  smallpot = cpot(keep, node_sizes(keep), bigpot.g, bigpot.h, bigpot.K);
else
  [h1, h2, K11, K12, K21, K22] = partition_matrix_vec(bigpot.h, bigpot.K, sum_over, keep, node_sizes);
  n = length(h1);
  K11inv = inv(K11);
  g = bigpot.g + 0.5*(n*log(2*pi) - log(det(K11)) + h1'*K11inv*h1);
  if length(h2) > 0 % ~isempty(keep) % we are are actually keeping something
    A = K21*K11inv;
    h = h2 - A*h1;
    K = K22 - A*K12;
  else
    h = [];
    K = [];
  end
  smallpot = cpot(keep, node_sizes(keep), g, h, K);
end
           
