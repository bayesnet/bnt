function [m,k,w] = vonMises_CPD_params_given_dps(CPD,domain,evidence)
%VONMISES_ Summary of this function goes here
%   Detailed explanation goes here
ps = domain(1:end-1);
dps = ps(CPD.dps);
if isempty(dps)
  m = CPD.mean;
  k = CPD.con;
  w = CPD.weights;
else
  odom = domain(~isemptycell(evidence(domain)));
  dops = myintersect(dps, odom);
  dpvals = cat(1, evidence{dops});
  if length(dops) == length(dps)
    dpsizes = CPD.sizes(CPD.dps);
    dpval = subv2ind(dpsizes, dpvals(:)');
    m = CPD.mean(:, dpval);
    k = CPD.con(:, :, dpval);
    w = CPD.weights(:,:,dpval);
  else
    map = find_equiv_posns(dops, dps);
    index = mk_multi_index(length(dps), map, dpvals);
    m = CPD.mean(:, index{:});
    k = CPD.con(:,:,index{:});
    w = CPD.weights(:,:,index{:});
  end
end

