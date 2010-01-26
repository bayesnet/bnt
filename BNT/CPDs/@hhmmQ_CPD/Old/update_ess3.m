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
% hor_counts(ndx{:}) = fmarginal(...)
% where e.g., ndx = {1, ':', 2} if Qps is hidden but we observe old_self=1, self=2.

% ndx{i,k,j}
if hidden_bitv(old_self)
  ndx{1} = ':';
else
  ndx{1} = evidence{old_self};
end
if hidden_bitv(Qps)
  ndx{2} = ':';
else
  ndx{2} = subv2ind(Qpsz, cat(1, evidence{Qps}));
end
if hidden_bitv(self)
  ndx{3} = ':';
else
  ndx{3} = evidence{self};
end

fmarg = add_ev_to_dmarginal(fmarginal, evidence, ns);
% marg(Qold(t-1), Fbelow(t-1), Fself(t-1), Qps(t), Qself(t))                            
hor_counts = zeros(Qsz, Qpsz, Qsz);
ver_counts = zeros(Qpsz, Qsz);
    
if ~isempty(CPD.Fbelow_ndx)
  if ~isempty(CPD.Fself_ndx) % general case
    fmarg.T = myreshape(fmarg.T, [Qsz 2 2 Qpsz Qsz]);
    marg_ndx = {ndx{1}, 2, 1, ndx{2}, ndx{3}};
    hor_counts(ndx{:}) = fmarg.T(marg_ndx{:});
    ver_counts(ndx{2:3}) = ... % sum over Fbelow and Qold=i
	sum(fmarg.T({ndx{1}, 1, 2, ndx{2}, ndx{3}}),1) + ..
	sum(fmarg.T({ndx{1}, 2, 2, ndx{2}, ndx{3}}),1);
  else % no F from self, hence no startprob
    fmarg.T = myreshape(fmarg.T, [Qsz 2 Qpsz Qsz]);
    hor_counts(ndx{:}) = fmarg.T({ndx{1}, 2, ndx{2}, ndx{3}});
  end
else % no F signal from below
  if ~isempty(CPD.Fself_ndx) % self F
    fmarg.T = myreshape(fmarg.T, [Qsz 2 Qpsz Qsz]);
    hor_counts(ndx{:}) = fmarg.T({ndx{1}, 1, ndx{2}, ndx{3}});
    ver_counts(ndx{2:3}) = ... % sum over Qold=i
	sum(fmarg.T({ndx{1}, 2, ndx{2}, ndx{3}}),1);
  else % no F from self
    error('An hhmmQ node without any F parents is just a tabular_CPD')
  end
end


CPD.sub_CPD_trans = update_ess_simple(CPD.sub_CPD_trans, hor_counts);

if ~isempty(CPD.sub_CPD_start)
  CPD.sub_CPD_start = update_ess_simple(CPD.sub_CPD_start, ver_counts);
end


