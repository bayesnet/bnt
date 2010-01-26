function CPD = bayes_update_params(CPD, self_ev, pev)
% UPDATE_PARAMS_COMPLETE Bayesian parameter updating given completely observed data (tabular)
% CPD = update_params_complete(CPD, self_ev, pev)
%
% self_ev(m) is the evidence on this node in case m.
% pev(i,m) is the evidence on the i'th parent in case m (if there are any parents).
% These can be arrays or cell arrays.
%
% We update the Dirichlet pseudo counts and set the CPT to the mean of the posterior.

if iscell(self_ev), usecell = 1; else usecell = 0; end

ncases = length(self_ev);
sz = CPD.sizes;
nparents = length(sz)-1;
assert(nparents == size(pev,1));

if ncases == 0 | ~adjustable_CPD(CPD)
  return;
elseif ncases == 1 % speedup the sequential learning case by avoiding normalization of the whole array
  if usecell
    x = cat(1, pev{:})';
    y = self_ev{1};
  else
    x = pev(:)';
    y = self_ev;
  end
  switch nparents
   case 0,
    CPD.dirichlet(y) = CPD.dirichlet(y)+1;
    CPD.CPT = CPD.dirichlet / sum(CPD.dirichlet);
   case 1,
    CPD.dirichlet(x(1), y) = CPD.dirichlet(x(1), y)+1;
    CPD.CPT(x(1), :) = CPD.dirichlet(x(1), :) ./ sum(CPD.dirichlet(x(1), :));
   case 2,
    CPD.dirichlet(x(1), x(2), y) = CPD.dirichlet(x(1), x(2), y)+1;
    CPD.CPT(x(1), x(2), :) = CPD.dirichlet(x(1), x(2), :) ./ sum(CPD.dirichlet(x(1), x(2), :));
   case 3,
    CPD.dirichlet(x(1), x(2), x(3), y) = CPD.dirichlet(x(1), x(2), x(3), y)+1;
    CPD.CPT(x(1), x(2), x(3), :) = CPD.dirichlet(x(1), x(2), x(3), :) ./ sum(CPD.dirichlet(x(1), x(2), x(3), :));
   otherwise,
    ind = subv2ind(sz, [x y]);
    CPD.dirichlet(ind) = CPD.dirichlet(ind) + 1;
    CPD.CPT = mk_stochastic(CPD.dirichlet);
  end
else  
  if usecell
    data = [cell2num(pev); cell2num(self_ev)]; 
  else
    data = [pev; self_ev];
  end
  counts = compute_counts(data, sz);
  CPD.dirichlet = CPD.dirichlet + counts;
  CPD.CPT = mk_stochastic(CPD.dirichlet);
end
