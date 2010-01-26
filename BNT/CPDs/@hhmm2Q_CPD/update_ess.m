function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, hidden_bitv)

marg = add_ev_to_dmarginal(fmarginal, evidence,  ns);

nps = length(CPD.dom_sz)-1; % num parents

if ~isempty(CPD.Fbelow_ndx)
  if ~isempty(CPD.Fself_ndx) % general case
    ndx = mk_multi_index(nps+1, [CPD.Fbelow_ndx CPD.Fself_ndx], [2 1]);
    CPD.trans_counts = CPD.trans_counts + squeeze(marg.T(ndx{:}));
    ndx = mk_multi_index(nps+1, [CPD.Fbelow_ndx CPD.Fself_ndx], [2 2]);
    CPD.start_counts = CPD.start_counts + squeeze(marg.T(ndx{:}));
  else % no F from self, hence no startprob (top level)
    ndx = mk_multi_index(nps+1, CPD.Fbelow_ndx, 2);
    CPD.trans_counts = CPD.trans_counts + squeeze(marg.T(ndx{:}));
  end
else % no F signal from below
  if ~isempty(CPD.Fself_ndx) % self F (bottom level)
    ndx = mk_multi_index(nps+1, CPD.Fself_ndx, 1);
    CPD.trans_counts = CPD.trans_counts + squeeze(marg.T(ndx{:}));
    ndx = mk_multi_index(nps+1, CPD.Fself_ndx, 2);
    CPD.start_counts = CPD.start_counts + squeeze(marg.T(ndx{:}));
  else % no F from self or below
    error('no F signal')
  end
end
