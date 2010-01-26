function T = convert_CPD_to_table_hidden_ps(CPD, child_obs)
% CONVERT_CPD_TO_TABLE_HIDDEN_PS Convert a discrete CPD to a table
% T = convert_CPD_to_table_hidden_ps(CPD, child_obs)
%
% This is like convert_to_table, except that we are guaranteed that
% none of the parents have evidence on them.
% child_obs may be an integer (1,2,...) or [].

CPT = CPD_to_CPT(CPD);
if isempty(child_obs)
  T = CPT(:);
else
  sz = dom_sizes(CPD);
  if length(sz)==1 % no parents
    T = CPT(child_obs);
  else
    CPT = reshape(CPT, prod(sz(1:end-1)), sz(end));
    T = CPT(:, child_obs);
  end
end
