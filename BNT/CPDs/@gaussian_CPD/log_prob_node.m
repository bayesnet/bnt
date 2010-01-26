function L = log_prob_node(CPD, self_ev, pev)
% LOG_PROB_NODE Compute prod_m log P(x(i,m)| x(pi_i,m), theta_i) for node i (gaussian)
% L = log_prob_node(CPD, self_ev, pev)
%
% self_ev(m) is the evidence on this node in case m.
% pev(i,m) is the evidence on the i'th parent in case m (if there are any parents).
% (These may also be cell arrays.)

if iscell(self_ev), usecell = 1; else usecell = 0; end

use_log = 1;
ncases = length(self_ev);
nparents = length(CPD.sizes)-1;
assert(ncases == size(pev, 2));

if ncases == 0
  L = 0;
  return;
end

L = 0;
for m=1:ncases
  if isempty(CPD.dps)
    i = 1;
  else
    if usecell
      dpvals = cat(1, pev{CPD.dps, m});
    else
      dpvals = pev(CPD.dps, m);
    end
    i = subv2ind(CPD.sizes(CPD.dps), dpvals(:)');
  end
  if usecell
    y = self_ev{m};
  else
    y = self_ev(m);
  end
  if length(CPD.cps) == 0 
    L = L + gaussian_prob(y, CPD.mean(:,i), CPD.cov(:,:,i), use_log);
  else
    if usecell
      x = cat(1, pev{CPD.cps, m});
    else
      x = pev(CPD.cps, m);
    end
    L = L + gaussian_prob(y, CPD.mean(:,i) + CPD.weights(:,:,i)*x, CPD.cov(:,:,i), use_log);
  end
end

