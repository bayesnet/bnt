function loglik = bethe_free_energy(engine, evidence)
% BETHE_FREE_ENERGY Compute Bethe free energy approximation to the log likelihood
% loglik = bethe_free_energy(engine, evidence)
%
% The Bethe free energy is given by an exact energy term and an approximate entropy term.
% Energy
%  E = -sum_f sum_i b(f,i) ln theta(f,i)
% where b(f,i) = approximate Pr(family f = i) 
% and theta(f,i) = Pr(f = i)
% Entropy
%  S = H1 - H2
%  H1 = sum_f sum_p H(b(f))
% where b(f) = belief on family f, H(.) = entropy
%  H2 = sum_n (q(n)-1) H(b(n))
% where q(n) = num. neighbors of n
%
% This function was written by Yair Weiss, 8/22/01.

hidden = find(isemptycell(evidence));
bnet = bnet_from_engine(engine);
N = length(bnet.dag);

add_ev = 1;
E=0;H1=0;H2=0;
loglik=0;
for n=1:N
  ps=parents(bnet.dag,n);
  if (length(ps)==0) % root node
    qi=length(children(bnet.dag,n))-1;
  else
    qi=length(children(bnet.dag,n));
  end
  bf = marginal_family(engine, n, add_ev);
  bf = bf.T(:);
  e = bnet.equiv_class(n);
  T = CPD_to_CPT(bnet.CPD{e});
  T = T(:);
  E = E-sum(log(T+(T==0)).*bf);

  if length(ps) > 0
    % root nodes don't count as fmailies
    H1 = H1+sum(log(bf+(bf==0)).*bf);
  end
  
  bi = marginal_nodes(engine, n, add_ev);
  bi = bi.T(:);
  H2 = H2+qi*sum(log(bi+(bi==0)).*bi);
end
loglik=E+H1-H2;
loglik=-loglik;
