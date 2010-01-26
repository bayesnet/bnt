function CPD = update_CPT(CPD)
% Compute the big CPT for an HHMM Q node (including F parents) given internal transprob and startprob
% function CPD = update_CPT(CPD)

Qsz = CPD.Qsz;
Qpsz = CPD.Qpsz;

if ~isempty(CPD.Fbelow_ndx)
  if ~isempty(CPD.Fself_ndx) % general case
    % Fb(t-1) Fself(t-1)  P(Q(t)=j| Q(t-1)=i, Qps(t)=k)
    % ------------------------------------------------------
    % 1        1         delta(i,j)
    % 2        1         transprob(i,k,j)
    % 1        2         impossible
    % 2        2         startprob(k,j)
    CPT = zeros(Qsz, 2, 2, Qpsz, Qsz);
    I = repmat(eye(Qsz), [1 1 Qpsz]); % i,j,k
    I = permute(I, [1 3 2]); % i,k,j
    CPT(:, 1, 1, :, :) = I;
    CPT(:, 2, 1, :, :) = CPD.transprob;
    CPT(:, 1, 2, :, :) = I;
    CPT(:, 2, 2, :, :) = repmat(reshape(CPD.startprob, [1 Qpsz Qsz]), [Qsz 1 1]); % replicate over i
  else % no F from self, hence no startprob
    % Fb(t-1) P(Q(t)=j| Q(t-1)=i, Qps(t)=k)
    % ------------------------------------------------------
    % 1       delta(i,j)
    % 2       transprob(i,k,j)
    
    nps = length(CPD.dom_sz)-1; % num parents
    CPT = 0*myones(CPD.dom_sz);
    %CPT = zeros(Qsz, 2, Qpsz, Qsz); % assumes CPT(Q(t-1), F(t-1), Qps, Q(t))
    % but a member of Qps may preceed Q(t-1) or F(t-1) in the ordering
    
    I = repmat(eye(Qsz), [1 1 Qpsz]); % i,j,k
    I = permute(I, [1 3 2]); % i,k,j

    % the following fails if there is a member of Qps with a lower
    % number than F
    %CPT(:, 1, :, :) = I;
    %CPT(:, 2, :, :) = CPD.transprob;

    ndx = mk_multi_index(nps+1, CPD.Fbelow_ndx, 1);
    CPT(ndx{:}) = I;
    ndx = mk_multi_index(nps+1, CPD.Fbelow_ndx, 2);
    CPT(ndx{:}) = CPD.transprob;
    keyboard
  end
else % no F signal from below
  if ~isempty(CPD.Fself_ndx)
    % Q(t-1), Fself(t-1), Qps, Q(t)
    
    % if condition start on previous concrete state (as in map learning),
    % CPT(:, 1, :, :, :) = CPD.transprob(Q(t-1), Qps, Q(t))
    % CPT(:, 2, :, :, :) = CPD.startprob(Q(t-1), Qps, Q(t))
    
    % Fself(t-1)  P(Q(t-1)=i, Qps(t)=k -> Q(t)=j)
    % ------------------------------------------------------
    % 1         transprob(i,k,j)
    % 2         startprob(k,j)
    CPT = zeros(Qsz, 2, Qpsz, Qsz);
    I = repmat(eye(Qsz), [1 1 Qpsz]); % i,j,k
    I = permute(I, [1 3 2]); % i,k,j
    CPT(:, 1, :, :) = CPD.transprob;
    if CPD.fullstartprob
      CPT(:, 2, :, :) = CPD.startprob;
    else
      CPT(:, 2, :, :) = repmat(reshape(CPD.startprob, [1 Qpsz Qsz]), [Qsz 1 1]); % replicate over i
    end
    else % no F from self
    error('An hhmmQ node without any F parents is just a tabular_CPD')
  end
end

CPD = set_fields(CPD, 'CPT', CPT);          
