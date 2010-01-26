function m = marginal_family(engine, i, t)
% MARGINAL_FAMILY Compute the marginal on i's family (inf_engine)
% m = marginal_family(engine, i, t)
%
% t defaults to 1.

if nargin < 3, t = 1; end

bnet = bnet_from_engine(engine);
if t==1
  m = marginal_nodes(engine, family(bnet.dag, i));
else
  ss = length(bnet.intra);
  fam = family(bnet.dag, i+ss);
  if any(fam<=ss)
    % i has a parent in the preceeding slice
    % Hence the lowest numbered slice containing the family is t-1
    m = marginal_nodes(engine, fam, t-1);
  else
    % The family all fits inside slice t
    % Hence shift the indexes back to slice 1
    m = marginal_nodes(engine, fam-ss, t);
  end
end     
