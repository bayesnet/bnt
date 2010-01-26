function L = log_marg_prob_node(CPD, self_ev, pev, usecell)
% LOG_MARG_PROB_NODE Compute sum_m log P(x(i,m)| x(pi_i,m)) for node i (tabular)
% L = log_marg_prob_node(CPD, self_ev, pev)
%
% This differs from log_prob_node because we integrate out the parameters.
% self_ev(m) is the evidence on this node in case m.
% pev(i,m) is the evidence on the i'th parent in case m (if there are any parents).
% (These may also be cell arrays.)

ncases = length(self_ev);
sz = CPD.sizes;
nparents = length(sz)-1;
assert(ncases == size(pev, 2)); 

if nargin < 4
  %usecell = 0;
  if iscell(self_ev)
    usecell = 1;
  else
    usecell = 0;
  end
end


if ncases==0
  L = 0;
  return;
elseif ncases==1  % speedup the sequential learning case
  CPT = CPD.CPT;
  % We assume the CPTs are already set to the mean of the posterior (due to bayes_update_params)
  if usecell
    x = cat(1, pev{:})';
    y = self_ev{1};
  else
    %x = pev(:)';
    x = pev;
    y = self_ev;
  end
  switch nparents
   case 0, p = CPT(y);
   case 1, p = CPT(x(1), y);
   case 2, p = CPT(x(1), x(2), y);
   case 3, p = CPT(x(1), x(2), x(3), y);
   otherwise,
    ind = subv2ind(sz, [x y]);
    p = CPT(ind);
  end
  L = log(p);
else
  % We ignore the CPTs here and assume the prior has not been changed
  
  % We arrange the data as in the following example.
  % Let there be 2 parents and 3 cases. Let p(i,m) be parent i in case m,
  % and y(m) be the child in case m. Then we create the data matrix
  % 
  % p(1,1) p(1,2) p(1,3)
  % p(2,1) p(2,2) p(2,3)
  % y(1)   y(2)   y(3)
  if usecell
    data = [cell2num(pev); cell2num(self_ev)]; 
  else
    data = [pev; self_ev];
  end
  %S = struct(CPD); fprintf('log marg prob node %d, ps\n', S.self); disp(S.parents)
  counts = compute_counts(data, sz);
  L = dirichlet_score_family(counts, CPD.dirichlet);
end


