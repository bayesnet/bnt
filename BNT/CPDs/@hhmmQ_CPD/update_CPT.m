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
    CPT(:, 2, 2, :, :) = repmat(reshape(CPD.startprob, [1 Qpsz Qsz]), ...
				[Qsz 1 1]); % replicate  over i 
  else % no F from self, hence no startprob
    % Fb(t-1) P(Q(t)=j| Q(t-1)=i, Qps(t)=k)
    % ------------------------------------------------------
    % 1       delta(i,j)
    % 2       transprob(i,k,j)
    
    nps = length(CPD.dom_sz)-1; % num parents
    CPT = 0*myones(CPD.dom_sz);
    %CPT = zeros(Qsz, 2, Qpsz, Qsz); % assumes CPT(Q(t-1), F(t-1), Qps, Q(t))
    % but a member of Qps may preceed Q(t-1) or F(t-1) in the ordering

    for k=1:CPD.Qpsz
      Qps_vals = ind2subv(CPD.Qpsizes, k);
      ndx = mk_multi_index(nps+1, [CPD.Fbelow_ndx CPD.Qps_ndx], [1 Qps_vals]);
      CPT(ndx{:}) = eye(Qsz); % CPT(:,2,k,:) or CPT(:,k,2,:) etc
    end
    ndx = mk_multi_index(nps+1, CPD.Fbelow_ndx, 2);
    CPT(ndx{:}) = CPD.transprob; % we assume transprob is in topo order
  end
else % no F signal from below
  if ~isempty(CPD.Fself_ndx)
    % Q(t-1), Fself(t-1), Qps, Q(t)
    
    % Fself(t-1)  P(Q(t-1)=i, Qps(t)=k -> Q(t)=j)
    % ------------------------------------------------------
    % 1         transprob(i,k,j)
    % 2         startprob(k,j)
    
    nps = length(CPD.dom_sz)-1; % num parents
    CPT = 0*myones(CPD.dom_sz);
    ndx = mk_multi_index(nps+1, CPD.Fself_ndx, 1);
    CPT(ndx{:}) = CPD.transprob;
    if CPD.fullstartprob
      ndx = mk_multi_index(nps+1, CPD.Fself_ndx, 2);
      CPT(ndx{:}) = CPD.startprob;
    else
      for i=1:CPD.Qsz
	ndx = mk_multi_index(nps+1, [CPD.Fself_ndx CPD.old_self_ndx], [2 i]);
	CPT(ndx{:}) = CPD.startprob;
      end
    end
  else % no F from self
    error('An hhmmQ node without any F parents is just a tabular_CPD')
  end
end

CPD = set_fields(CPD, 'CPT', CPT);          
