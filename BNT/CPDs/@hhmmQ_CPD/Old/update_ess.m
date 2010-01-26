function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, hidden_bitv)
% UPDATE_ESS Update the Expected Sufficient Statistics of a hhmm Q node.
% function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, idden_bitv)

% Figure out the node numbers associated with each parent
% e.g., D=4, d=3, Qps = all Qs above, so dom = [Q3(t-1) F4(t-1) F3(t-1) Q1(t) Q2(t) Q3(t)].
% so self = Q3(t), old_self = Q3(t-1), CPD.Qps = [1 2], Qps = [Q1(t) Q2(t)]
dom = fmarginal.domain;
self = dom(end);
old_self = dom(1);
Qps = dom(length(dom)-length(CPD.Qps):end-1);

Qsz = CPD.Qsizes(CPD.d);
Qpsz = prod(CPD.Qsizes(CPD.Qps));

% If some of the Q nodes are observed (which happens during supervised training)
% the counts will only be non-zero in positions
% consistent with the evidence. We put the computed marginal responsibilities
% into the appropriate slots of the big counts array.
% (Recall that observed discrete nodes only have a single effective value.)
% (A more general, but much slower, way is to call add_evidence_to_dmarginal.)
% We assume the F nodes are never observed.

obs_self = ~hidden_bitv(self);
obs_Qps = (~isempty(Qps)) & (~any(hidden_bitv(Qps))); % we assume that all or none of the Q parents are observed

if obs_self
  self_val = evidence{self};
  oldself_val = evidence{old_self};
end

if obs_Qps
  Qps_val = subv2ind(Qpsz, cat(1, evidence{Qps}));
  if Qps_val == 0
    keyboard
  end
end

if CPD.d==1 % no Qps from above
  if ~CPD.F1toQ1 % no F from self
    % marg(Q1(t-1), F2(t-1), Q1(t))                            
    % F2(t-1) P(Q1(t)=j | Q1(t-1)=i)
    % 1       delta(i,j)
    % 2       transprob(i,j)
    if obs_self
      hor_counts = zeros(Qsz, Qsz);
      hor_counts(oldself_val, self_val) = fmarginal.T(2);
    else
      marg = reshape(fmarginal.T, [Qsz 2 Qsz]);
      hor_counts = squeeze(marg(:,2,:));
    end
  else
    % marg(Q1(t-1), F2(t-1), F1(t-1), Q1(t))                            
    % F2(t-1) F1(t-1)  P(Qd(t)=j| Qd(t-1)=i)
    % ------------------------------------------------------
    % 1        1         delta(i,j)
    % 2        1         transprob(i,j)
    % 1        2         impossible
    % 2        2         startprob(j)
    if obs_self
      marg = myreshape(fmarginal.T, [1 2 2 1]);
      hor_counts = zeros(Qsz, Qsz);
      hor_counts(oldself_val, self_val) = marg(1,2,1,1);
      ver_counts = zeros(Qsz, 1);
      %ver_counts(self_val) = marg(1,2,2,1);
      ver_counts(self_val) = marg(1,2,2,1) + marg(1,1,2,1);
    else
      marg = reshape(fmarginal.T, [Qsz 2 2 Qsz]);
      hor_counts = squeeze(marg(:,2,1,:));
      %ver_counts = squeeze(sum(marg(:,2,2,:),1)); % sum over i
      ver_counts = squeeze(sum(marg(:,2,2,:),1)) + squeeze(sum(marg(:,1,2,:),1)); % sum i,b
    end
  end % F1toQ1
else % d ~= 1
  if CPD.d < CPD.D % general case
    % marg(Qd(t-1), Fd+1(t-1), Fd(t-1), Qps(t), Qd(t))                            
    % Fd+1(t-1) Fd(t-1)  P(Qd(t)=j| Qd(t-1)=i, Qps(t)=k)
    % ------------------------------------------------------
    % 1        1         delta(i,j)
    % 2        1         transprob(i,k,j)
    % 1        2         impossible
    % 2        2         startprob(k,j)
    if obs_Qps & obs_self
      marg = myreshape(fmarginal.T, [1 2 2 1 1]);
      k = 1;
      hor_counts = zeros(Qsz, Qpsz, Qsz);
      hor_counts(oldself_val, Qps_val, self_val) = marg(1, 2,1, k,1);
      ver_counts = zeros(Qpsz, Qsz);
      %ver_counts(Qps_val, self_val) = marg(1, 2,2, k,1);
      ver_counts(Qps_val, self_val) = marg(1, 2,2, k,1) + marg(1, 1,2, k,1);
    elseif obs_Qps & ~obs_self
      marg = myreshape(fmarginal.T, [Qsz 2 2 1 Qsz]);
      k = 1;
      hor_counts = zeros(Qsz, Qpsz, Qsz);
      hor_counts(:, Qps_val, :) = marg(:, 2,1, k,:);
      ver_counts = zeros(Qpsz, Qsz);
      %ver_counts(Qps_val, :) = sum(marg(:, 2,2, k,:), 1);
      ver_counts(Qps_val, :) = sum(marg(:, 2,2, k,:), 1) + sum(marg(:, 1,2, k,:), 1);
    elseif ~obs_Qps & obs_self
      error('not yet implemented')
    else % everything is hidden
      marg = reshape(fmarginal.T, [Qsz 2 2 Qpsz Qsz]);
      hor_counts = squeeze(marg(:,2,1,:,:)); % i,k,j
      %ver_counts = squeeze(sum(marg(:,2,2,:,:),1)); % sum over i
      ver_counts = squeeze(sum(marg(:,2,2,:,:),1)) + squeeze(sum(marg(:,1,2,:,:),1)); % sum over i,b
    end
  else % d == D, so no F from below
    % marg(QD(t-1), FD(t-1), Qps(t), QD(t))                            
    % FD(t-1) P(QD(t)=j | QD(t-1)=i, Qps(t)=k)
    % 1      transprob(i,k,j) 
    % 2      startprob(k,j)
    if obs_Qps & obs_self
      marg = myreshape(fmarginal.T, [1 2 1 1]);
      k = 1;
      hor_counts = zeros(Qsz, Qpsz, Qsz);
      hor_counts(oldself_val, Qps_val, self_val) = marg(1, 1, k,1);
      ver_counts = zeros(Qpsz, Qsz);
      ver_counts(Qps_val, self_val) = marg(1, 2, k,1);
    elseif obs_Qps & ~obs_self
      marg = myreshape(fmarginal.T, [Qsz 2 1 Qsz]);
      k = 1;
      hor_counts = zeros(Qsz, Qpsz, Qsz);
      hor_counts(:, Qps_val, :) = marg(:, 1, k,:);
      ver_counts = zeros(Qpsz, Qsz);
      ver_counts(Qps_val, :) = sum(marg(:, 2, k, :), 1);
    elseif ~obs_Qps & obs_self
      error('not yet implemented')
    else % everything is hidden
      marg = reshape(fmarginal.T, [Qsz 2 Qpsz Qsz]);
      hor_counts = squeeze(marg(:,1,:,:));
      ver_counts = squeeze(sum(marg(:,2,:,:),1)); % sum over i
    end
  end
end

CPD.sub_CPD_trans = update_ess_simple(CPD.sub_CPD_trans, hor_counts);

if ~isempty(CPD.sub_CPD_start)
  CPD.sub_CPD_start = update_ess_simple(CPD.sub_CPD_start, ver_counts);
end

