function CPT = CPD_to_CPT(CPD)
% Compute the big CPT for an HHMM Q node (including F parents)
% by combining internal transprob and startprob
% function CPT = CPD_to_CPT(CPD)

Qsz = CPD.Qsz;

if ~isempty(CPD.Fbelow_ndx)
  if ~isempty(CPD.Fself_ndx) % general case
    error('not implemented')
  else % no F from self, hence no startprob (top level)
    nps = length(CPD.dom_sz)-1; % num parents
    CPT = 0*myones(CPD.dom_sz);
    % when Fself=1, the CPT(i,j) = delta(i,j) for all k
    for k=1:prod(CPD.Qpsizes)
      Qps_vals = ind2subv(CPD.Qpsizes, k);
      ndx = mk_multi_index(nps+1, [CPD.Fbelow_ndx CPD.Qps_ndx], [1 Qps_vals]);
      CPT(ndx{:}) = eye(Qsz); % CPT(:,2,k,:) or CPT(:,k,2,:) etc
    end
    ndx = mk_multi_index(nps+1, CPD.Fbelow_ndx, 2);
    CPT(ndx{:}) = CPD.transprob; % we assume transprob is in topo order
  end
else % no F signal from below
  if ~isempty(CPD.Fself_ndx) % bottom level
    nps = length(CPD.dom_sz)-1; % num parents
    CPT = 0*myones(CPD.dom_sz);
    ndx = mk_multi_index(nps+1, CPD.Fself_ndx, 1);
    CPT(ndx{:}) = CPD.transprob;
    ndx = mk_multi_index(nps+1, CPD.Fself_ndx, 2);
    CPT(ndx{:}) = CPD.startprob;
  else % no F from self
    error('An hhmmQ node without any F parents is just a tabular_CPD')
  end
end

