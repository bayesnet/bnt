function strides = compute_strides(bnet)
% COMPUTE_STRIDES For each CPT and each variable in that CPT,
% returns the stride of that variable.  So in future, we can
% quickly extract a slice of the CPT.
%
% The return value is a 2d array, where strides(i,j) contains the
% stride of the jth variable in the ith CPT.  Cell arrays would
% have saved space but they are slower.
% 

num_cpts = size(bnet.CPD, 2);
max_cpt_dim = 1 + max(sum(bnet.dag));
strides = zeros(num_cpts, max_cpt_dim);

for i = 1:num_cpts
  c = CPT(bnet, i);
  siz = size(CPT(bnet, i));
  
  % Deal with the special case of a 1-d array separately
  if siz(2) == 1
    dim = 1;
  else
    dim = size(siz, 2);
  end

  strides(i, 1:dim ) = [1 cumprod(siz(1:dim-1))];
end
