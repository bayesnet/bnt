function T = convert_to_table(CPD, domain, local_ev, obs_bitv)
% CONVERT_TO_TABLE Convert a discrete CPD to a table
% function T = convert_to_table(CPD, domain, local_ev, obs_bitv)
%
% We convert the CPD to a CPT, and then lookup the evidence on the discrete parents.
% The resulting table can easily be converted to a potential.


CPT = CPD_to_CPT(CPD);
obs_child_only = ~any(obs_bitv(1:end-1)) & obs_bitv(end);

if obs_child_only
  sz = size(CPT);
  CPT = reshape(CPT, prod(sz(1:end-1)), sz(end));
  o = local_ev{end};
  T = CPT(:, o);
else
  odom = domain(obs_bitv);  
  vals = cat(1, local_ev{find(obs_bitv)}); % undo cell array
  map = find_equiv_posns(odom, domain);
  index = mk_multi_index(length(domain), map, vals);
  T = CPT(index{:});
end
