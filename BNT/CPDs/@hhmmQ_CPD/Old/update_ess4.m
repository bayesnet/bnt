function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, hidden_bitv)
% UPDATE_ESS Update the Expected Sufficient Statistics of a hhmm Q node.
% function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, idden_bitv)
%
% we assume if one of the Qps is observed, all of them are 
% We assume the F nodes are already hidden 

% Figure out the node numbers associated with each parent
dom = fmarginal.domain;
self = dom(CPD.self_ndx);
old_self = dom(CPD.old_self_ndx);
%Fself = dom(CPD.Fself_ndx);
%Fbelow = dom(CPD.Fbelow_ndx);
Qps = dom(CPD.Qps_ndx);

Qsz = CPD.Qsz;
Qpsz = CPD.Qpsz;


% hor_counts(old_self, Qps, self),
% fmarginal(old_self, Fbelow, Fself, Qps, self)
% hor_counts(i,k,j) = fmarginal(i,2,1,k,j) % below has finished, self has not
% ver_counts(i,k,j) = fmarginal(i,2,2,k,j) % below has finished, and so has self (reset)
% Since any of i,j,k may be observed, we write
% hor_counts(i_counts_ndx, kndx, jndx) = fmarginal(i_fmarg_ndx...)
% where i_fmarg_ndx = 1 and i_counts_ndx = i if old_self is observed to have value i,
% i_fmarg_ndx = 1:Qsz and i_counts_ndx = 1:Qsz if old_self is hidden, etc.


if hidden_bitv(old_self)
  i_counts_ndx = 1:Qsz;
  i_fmarg_ndx = 1:Qsz;
  eff_oldQsz = Qsz;
else
  i_counts_ndx = evidence{old_self};
  i_fmarg_ndx = 1;
  eff_oldQsz = 1;
end

if all(hidden_bitv(Qps)) % we assume all are hidden or all are observed
  k_counts_ndx = 1:Qpsz;
  k_fmarg_ndx = 1:Qpsz;
  eff_Qpsz = Qpsz;
else
  k_counts_ndx = subv2ind(Qpsz, cat(1, evidence{Qps}));
  k_fmarg_ndx = 1;
  eff_Qpsz = 1;
end

if hidden_bitv(self)
  j_counts_ndx = 1:Qsz;
  j_fmarg_ndx = 1:Qsz;
  eff_Qsz = Qsz;
else
  j_counts_ndx = evidence{self};
  j_fmarg_ndx = 1;
  eff_Qsz = 1;
end

hor_counts = zeros(Qsz, Qpsz, Qsz);
ver_counts = zeros(Qpsz, Qsz);
    
if ~isempty(CPD.Fbelow_ndx)
  if ~isempty(CPD.Fself_ndx) % general case
    fmarg.T = myreshape(fmarg.T, [eff_oldQsz 2 2 eff_Qpsz eff_Qsz]);
    hor_counts(i_counts_ndx, k_counts_ndx, j_counts_ndx) = ...
	fmarg.T(:, i_fmarg_ndx, 2, 1, k_fmarg_ndx, j_fmarg_ndx);
    ver_counts(k_counts_ndx, j_counts_ndx) = ... % sum over Fbelow and Qold
	sum(fmarg.T(:, 1, 2, k_fmarg_ndx, j_fmarg_ndx), 1) + ...
	sum(fmarg.T(:, 2, 2, k_fmarg_ndx, j_fmarg_ndx), 1); 
  else % no F from self, hence no startprob
    fmarg.T = myreshape(fmarg.T, [eff_oldQsz 2 eff_Qpsz eff_Qsz]);
    hor_counts(i_counts_ndx, k_counts_ndx, j_counts_ndx) = ...
	fmarg.T(i_fmarg_ndx, 2, k_fmarg_ndx, j_fmarg_ndx);
  end
else % no F signal from below
  if ~isempty(CPD.Fself_ndx) % self F
    fmarg.T = myreshape(fmarg.T, [eff_oldQsz 2 eff_Qpsz eff_Qsz]);
    hor_counts(i_counts_ndx, k_counts_ndx, j_counts_ndx) = ...
	fmarg.T(i_fmarg_ndx, 1, k_fmarg_ndx, j_fmarg_ndx);
    ver_counts(k_counts_ndx, j_counts_ndx) = ... % sum over Qold
	sum(fmarg.T(:, 2, k_fmarg_ndx, j_fmarg_ndx), 1);
  else % no F from self
    error('An hhmmQ node without any F parents is just a tabular_CPD')
  end
end


CPD.sub_CPD_trans = update_ess_simple(CPD.sub_CPD_trans, hor_counts);

if ~isempty(CPD.sub_CPD_start)
  CPD.sub_CPD_start = update_ess_simple(CPD.sub_CPD_start, ver_counts);
end


