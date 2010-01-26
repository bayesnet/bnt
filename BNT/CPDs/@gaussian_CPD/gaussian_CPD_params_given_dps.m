function [m, C, W] = gaussian_CPD_params_given_dps(CPD, domain, evidence)
% GAUSSIAN_CPD_PARAMS_GIVEN_EV_ON_DPS Extract parameters given evidence on all discrete parents
% function [m, C, W] = gaussian_CPD_params_given_ev_on_dps(CPD, domain, evidence)

ps = domain(1:end-1);
dps = ps(CPD.dps);
if isempty(dps)
  m = CPD.mean;
  C = CPD.cov;
  W = CPD.weights;
else
  odom = domain(~isemptycell(evidence(domain)));
  dops = myintersect(dps, odom);
  dpvals = cat(1, evidence{dops});
  if length(dops) == length(dps)
    dpsizes = CPD.sizes(CPD.dps);
    dpval = subv2ind(dpsizes, dpvals(:)');
    m = CPD.mean(:, dpval);
    C = CPD.cov(:, :, dpval);
    W = CPD.weights(:, :, dpval);
  else
    map = find_equiv_posns(dops, dps);
    index = mk_multi_index(length(dps), map, dpvals);
    m = CPD.mean(:, index{:});
    C = CPD.cov(:, :, index{:});
    W = CPD.weights(:, :, index{:});
  end
end
