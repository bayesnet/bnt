function [D,d]= compute_bnet_nparams(bnet)
% [D,d] = compute_bnet_nparams(bnet)
%
% D is the dimension of the network
% d is the vector containing the number of parameters for all nodes

N = length(bnet.dag);
d = zeros(1,N);
for i=1:N
   a = struct(bnet.CPD{i});
   d(i) = a.nparams;
end
D = sum(d);
