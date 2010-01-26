function [engine, loglik] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (cond_gauss)
% [engine, loglik] = enter_evidence(engine, evidence, ...)
%
% evidence{i} = [] if if X(i) is hidden, and otherwise contains its observed value (scalar or column vector)

bnet = bnet_from_engine(engine);
ns = bnet.node_sizes(:);
observed = ~isemptycell(evidence);
onodes = find(observed);
hnodes = find(isemptycell(evidence));
engine.evidence = evidence;

% check there are no C->D links where C is hidden
pot_type = determine_pot_type(bnet, onodes);

dhid = myintersect(hnodes, bnet.dnodes);
S = prod(ns(dhid));
T = zeros(S,1);

N = length(bnet.dag);
mu = cell(1,N);
Sigma = cell(1,N); 
cobs = myintersect(bnet.cnodes, onodes);
chid = myintersect(bnet.cnodes, hnodes);
ens = ns;
ens(cobs) = 0;
for j=chid(:)'
  mu{j} = zeros(ens(j), S);
  Sigma{j} = zeros(ens(j), ens(j), S);
end
 
for i=1:S
  dvals = ind2subv(ns(dhid), i);
  evidence(dhid) = num2cell(dvals);
  [sub_engine, loglik] = enter_evidence(engine.sub_engine, evidence);
  for j=chid(:)'
    m = marginal_nodes(sub_engine, j);
    mu{j}(:,i) = m.mu;
    Sigma{j}(:,:,i) = m.Sigma;
  end
  T(i) = exp(loglik);
end

[T, lik] = normalise(T);
loglik = log(lik);

engine.T = T;
engine.mu = mu;
engine.Sigma = Sigma;

dnodes = bnet.dnodes;
dobs = myintersect(dnodes, onodes);
ens(dobs) = 1;
engine.joint_dmarginal = dpot(dnodes, ens(dnodes), myreshape(engine.T, ens(dnodes)));

engine.onodes = onodes;
