function T = convert_CPD_to_table_hidden_ps(CPD, self_val)
% CONVERT_CPD_TO_TABLE_HIDDEN_PS Convert a Gaussian CPD to a table
% function T = convert_CPD_to_table_hidden_ps(CPD, self_val)
%
% self_val must be a non-empty vector.
% All the parents are hidden.
%
% This is used by misc/convert_dbn_CPDs_to_tables

m = CPD.mean;
C = CPD.cov;
W = CPD.weights;

[ssz dpsize] = size(m);

T = zeros(dpsize, 1);
for i=1:dpsize
  T(i) = gaussian_prob(self_val, m(:,i), C(:,:,i));
end

