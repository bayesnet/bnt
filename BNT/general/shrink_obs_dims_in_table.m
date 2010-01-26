function T2 = shrink_obs_dims_in_table(T1, dom, evidence)
% SHRINK_OBS_DIMS_IN_TABLE Set observed dimensions to size 1
% T2 = shrink_obs_dims_in_table(T1, dom, evidence)
%
% If 'T1' contains observed nodes, it will have 0s in the positions that are
% inconsistent with the evidence. We now remove these 0s and set the corresponding dimensions to
% size 1, to be consistent with the way most inference engines handle evidence, which is to
% shrink observed nodes before doing inference.

% This is used by pearl and enumerative inf. engines.

odom = dom(~isemptycell(evidence(dom)));
vals = cat(1,evidence{odom});
ndx = mk_multi_index(length(dom), find_equiv_posns(odom, dom), vals(:));
T2 = T1(ndx{:});
