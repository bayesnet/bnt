function [transprob, termprob] = remove_hhmm_end_state(A)
% REMOVE_END_STATE Infer transition and termination probabilities from automaton with an end state
% [transprob, termprob] = remove_end_state(A)
% A(i,k,j) = Pr( i->j | Qps=k), where i in 1:Q, j in 1:(Q+1), and Q+1 is the end state

if ndims(A)==2 % top level
  Q = size(A,1);
  transprob = A(:,1:Q);
  termprob = A(:,Q+1)';
  
  % rescale
  for i=1:Q
    for j=1:Q
      denom = (1-termprob(i));
      denom = denom + (denom==0)*eps;
      transprob(i,j) = transprob(i,j) / denom;
    end
  end    
else
  Q = size(A,1);
  Qk = size(A,2);
  transprob = A(:, :, 1:Q);
  termprob = A(:,:,Q+1)';

  % rescale
  for k=1:Qk
    for i=1:Q
      for j=1:Q
	denom = (1-termprob(k,i));
	denom = denom + (denom==0)*eps;
	transprob(i,k,j) = transprob(i,k,j) / denom;
      end
    end    
  end
  
end

