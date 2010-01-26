function T = convert_to_table(CPD, domain, evidence)
% CONVERT_TO_TABLE Convert a discrete CPD to a table
% T = convert_to_table(CPD, domain, evidence)
%
% We convert the CPD to a CPT, and then lookup the evidence on the discrete parents.
% The resulting table can easily be converted to a potential.

CPT = CPD_to_CPT(CPD);
odom = domain(~isemptycell(evidence(domain)));
vals = cat(1, evidence{odom});
map = find_equiv_posns(odom, domain);
index = mk_multi_index(length(domain), map, vals);
T = CPT(index{:});
