function dag = sample_dag(P)
% SAMPLE_DAG Create a random directed acyclic graph with edge probabilities P(i,j)
% dag = sample_dag(P)
%
% This uses rejection sampling to reject graphs with directed cycles.

done = 0;
directed = 1;
iter = 1;
while ~done
  dag = binornd(1, P); % each edge is an indep Bernoulli (0/1) random variable
  dag = setdiag(dag, 0);
  done = acyclic(dag, directed);
  iter = iter + 1
end
