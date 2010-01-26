function [transprob, termprob] = remove_hhmm_end_state(A)
% REMOVE_END_STATE Infer transition and termination probabilities from automaton with an end state
% [transprob, termprob] = remove_end_state(A)
%
% A(i,k,j) = Pr( i->j | Qps=k), where i in 1:Q, j in 1:(Q+1), and Q+1 is the end state
% This implements the equation in footnote 3 of my NIPS 01 paper,
% transprob(i,k,j) = \tilde{A}_k(i,j)
% termprob(k,j) = \tau_k(j)
%
% For the top level, the k index is missing.

Q = size(A,1);
toplevel = (ndims(A)==2);
if toplevel
  Qk = 1;
  A = reshape(A, [Q 1 Q+1]);
else
  Qk = size(A, 2);
end

transprob = A(:, :, 1:Q);
term = A(:,:,Q+1)'; % term(k,j) = P(Qj -> end | k)
termprob = term;
%termprob = zeros(Qk, Q, 2);
%termprob(:,:,2) = term;
%termprob(:,:,1) = 1-term;

for k=1:Qk
  for i=1:Q
    for j=1:Q
      denom = (1-termprob(k,i));
      denom = denom + (denom==0)*eps;
      transprob(i,k,j) = transprob(i,k,j) / denom;
    end
  end    
end

if toplevel
  termprob = squeeze(termprob);
  transprob = squeeze(transprob);
end
