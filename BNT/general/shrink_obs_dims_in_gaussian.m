function marg2 = shrink_obs_dims_in_gaussian(marg1, dom, evidence, ns)
% SHRINK_OBS_DIMS_IN_GAUSSIAN Remove observed dimensions from mu/Sigma
% function marg2 = shrink_obs_dims_in_gaussian(marg1, dom, evidence, ns)

% This is used by loopy

hdom = dom(isemptycell(evidence(dom)));
ndx = find_equiv_posns(hdom, dom);
b = block(ndx, ns(dom));
marg2.mu = marg1.mu(b);
marg2.Sigma = marg1.Sigma(b,b);
marg2.domain = marg1.domain;
