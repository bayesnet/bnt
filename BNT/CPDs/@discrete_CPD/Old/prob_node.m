function [P, p] = prob_node(CPD, self_ev, pev)
% PROB_NODE Compute prod_m P(x(i,m)| x(pi_i,m), theta_i) for node i (discrete)
% [P, p] = prob_node(CPD, self_ev, pev)
%
% self_ev(m) is the evidence on this node in case m.
% pev(i,m) is the evidence on the i'th parent in case m (if there are any parents).
% (These may also be cell arrays.)
%
% p(m) = P(x(i,m)| x(pi_i,m), theta_i) 
% P = prod p(m)

if iscell(self_ev), usecell = 1; else usecell = 0; end

ncases = length(self_ev);
sz = dom_sizes(CPD);

nparents = length(sz)-1;
if nparents == 0
  assert(isempty(pev));
else
  assert(isequal(size(pev), [nparents ncases]));
end

n = length(sz);
dom = 1:n;
p = zeros(1, ncases);
if nparents == 0
  for m=1:ncases
    if usecell
      evidence = {self_ev{m}};
    else
      evidence = num2cell(self_ev(m));
    end
    T = convert_to_table(CPD, dom, evidence);
    p(m) = T;
  end
else
  for m=1:ncases
    if usecell
      evidence = cell(1,n);
      evidence(1:n-1) = pev(:,m);
      evidence(n) = self_ev(m);
    else
      evidence = num2cell([pev(:,m)', self_ev(m)]);
    end
    T = convert_to_table(CPD, dom, evidence);
    p(m) = T;
  end
end
P = prod(p);

