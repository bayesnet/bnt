function G = mk_adjmat_chain(T)
% MK_DAG_CHAIN Make adjacency matrix for bi-directional Markov chain of  T nodes
% function G = mk_dag_chain(T)
%
% G(t,t+1) = 1 for all t<T

G = diag(ones(1,T-1),1) + diag(ones(1,T-1),-1);

