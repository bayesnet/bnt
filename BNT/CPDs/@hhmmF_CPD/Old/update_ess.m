function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, hidden_bitv)
% UPDATE_ESS Update the Expected Sufficient Statistics of a hhmmF node.
% function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, hidden_bitv)

% Figure out the node numbers associated with each parent
% so we extract evidence from the right place
dom = fmarginal.domain; % Q(1) .. Q(d) F(d+1) F(d)
Qps = fmarginal.domain(1:end-2);
Q = Qps(end);
Qps = Qps(1:end-1);

Qsz = CPD.Qsizes(CPD.Q);
Qpsz = prod(CPD.Qsizes(CPD.Qps)); % may be 1

% We assume the F node are always hidden, but allow some of the Q nodes
% to be observed. We do case analysis for speed.
%We only extract prob from fmarginal.T when F(d+1)=2 i.e., model below has finished.
% wrong -> % We sum over the possibilities that F(d+1) = 1 or 2

obs_self = ~hidden_bitv(Q);
if obs_self
  self_val = evidence{Q};
end

if isempty(Qps) % independent of parent context
  counts = zeros(Qsz, 2);
  %fmarginal.T(Q(d), F(d+1), F(d))
  if obs_self
    marg = myreshape(fmarginal.T, [1 2 2]);
    counts(self_val,:) = marg(1,2,:);
    %counts(self_val,:) = marg(1,1,:) + marg(1,2,:);
  else
    marg = myreshape(fmarginal.T, [Qsz 2 2]);
    counts = squeeze(marg(:,2,:));
    %counts = squeeze(marg(:,2,:)) + squeeze(marg(:,1,:));
  end
else
  counts = zeros(Qpsz, Qsz, 2);
  %fmarginal.T(Q(1:d-1), Q(d), F(d+1), F(d))
  obs_Qps = ~any(hidden_bitv(Qps));  % we assume that all or none of the Q  parents are observed
  if obs_Qps
    Qps_val = subv2ind(Qpsz, cat(1, evidence{Qps}));
  end
  if obs_self & obs_Qps
    marg = myreshape(fmarginal.T, [1 1 2 2]);
    counts(Qps_val, self_val, :) = squeeze(marg(1,1,2,:));
    %counts(Qps_val, self_val, :) = squeeze(marg(1,1,2,:)) + squeeze(marg(1,1,1,:));
  elseif ~obs_self & obs_Qps
    marg = myreshape(fmarginal.T, [1 Qsz 2 2]);
    counts(Qps_val, :, :) = squeeze(marg(1,:,2,:));
    %counts(Qps_val, :, :) = squeeze(marg(1,:,2,:)) + squeeze(marg(1,:,1,:));
  elseif obs_self & ~obs_Qps
    error('not yet implemented')
  else
    marg = myreshape(fmarginal.T, [Qpsz Qsz 2 2]);
    counts(:, :, :) = squeeze(marg(:,:,2,:));
    %counts(:, :, :) = squeeze(marg(:,:,2,:)) + squeeze(marg(:,:,1,:));
  end    
end

CPD.sub_CPD_term = update_ess_simple(CPD.sub_CPD_term, counts);
