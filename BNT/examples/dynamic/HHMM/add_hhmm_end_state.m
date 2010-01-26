function A = add_hhmm_end_state(transprob, termprob)
% ADD_HMM_END_STATE Combine trans and term probs into transmat for automaton with an end state
% function A = add_hhmm_end_state(transprob, termprob)
%
% A(i,k,j) = Pr( i->j | Qps=k), where i in 1:Q, j in 1:(Q+1), and Q+1 is the end state
% This implements the equation in sec 4.6 of my tech report, where
% transprob(i,k,j) = \tilde{A}_k(i,j), termprob(k,j) = \tau_k(j)
%
% For the top level, the k index is missing.

Q = size(transprob,1);
toplevel = (ndims(transprob)==2);
if toplevel
  Qk = 1;
  transprob = reshape(transprob, [Q 1 Q]);
  termprob = reshape(termprob, [1 Q]);
else
  Qk = size(transprob, 2);
end

A = zeros(Q, Qk, Q+1);
A(:,:,Q+1) = termprob';

for k=1:Qk
  for i=1:Q
    for j=1:Q
      A(i,k,j) = transprob(i,k,j) * (1-termprob(k,i));
    end
  end    
end

if toplevel
  A = squeeze(A);
end
