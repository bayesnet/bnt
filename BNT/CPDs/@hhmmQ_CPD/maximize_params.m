function CPD = maximize_params(CPD, temp)
% MAXIMIZE_PARAMS Set the params of a hhmmQ node to their ML/MAP values.
% CPD = maximize_params(CPD, temperature)

Qsz = CPD.Qsz;
Qpsz = CPD.Qpsz;

if ~isempty(CPD.sub_CPD_start)
  CPD.sub_CPD_start = maximize_params(CPD.sub_CPD_start, temp);
  S = struct(CPD.sub_CPD_start);
  CPD.startprob = myreshape(S.CPT, [Qpsz Qsz]);
  %CPD.startprob = S.CPT;
end

if 1
  % If we are in a state that can only go the end state,
  % we will never see a transition to another (non-end) state,
  % so counts(i,k,j)=0 (and termprob(k,i)=1).
  % We set counts(i,k,i)=1 in this case.
  % This will cause remove_hhmm_end_state to return a
  % stochastic matrix, but otherwise has no effect on EM.
  counts = get_field(CPD.sub_CPD_trans, 'counts');
  counts = reshape(counts, [Qsz Qpsz Qsz]);
  for k=1:Qpsz
    for i=1:Qsz
      if sum(counts(i,k,:))==0 % never witnessed a transition out of i
	counts(i,k,i)=1; % add self loop 
	%fprintf('CPDQ d=%d i=%d k=%d\n', CPD.d, i, k);
      end
    end
  end
  CPD.sub_CPD_trans = set_fields(CPD.sub_CPD_trans, 'counts', counts(:)); 
end
 
CPD.sub_CPD_trans = maximize_params(CPD.sub_CPD_trans, temp);
S = struct(CPD.sub_CPD_trans);
%CPD.transprob = S.CPT;
CPD.transprob = myreshape(S.CPT, [Qsz Qpsz Qsz]);

CPD = update_CPT(CPD);
