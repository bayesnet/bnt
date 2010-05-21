function [n_dags, dag_list] = Markov_equivalent_dags(dag)

%
% [n_dags, dag_list] = Markov_equivalent_dags(dag) 
%
% generates a cell array of all Markov equivalent DAGs
% corresponding to the input DAG.
%
% YOU NEED TO HAVE THE STRUCTURE LEARNING PACKAGE IN PLACE TO USE THIS FUNCTION!
%
% Input:  DAG (in standard format, i.e. dag(a,b)=1 if and only if a->b)
%
% Output: Number of DAGs generated and 
%         Cell array of all Markov-equivalent DAGs (in same format as input)
% 
% Sample Use:
%
%   Example 1:  
%     % Find all DAGs equivalent to DAG of Asia Network
%     BN = mk_asia_bnet();
%     dag = BN.dag;
%     [n_dags, dag_list] = Markov_equivalent_dags(dag);
%     n_dags         % Answer should be: 6
%     dag_list{1}  % displays the first DAG, etc.
%
%   Example 2:  
%     % Find all DAGs equivalent to random DAG
%     dag = mk_rnd_dag(4);
%     [n_dags, dag_list] = Markov_equivalent_dags(dag);
%
% Imme Ebert-Uphoff (ebert@tree.com), 2007
%

  % find completed PDAG corresponding to DAG
  cpdag = dag_to_cpdag(dag);

  % convert to our notation, i.e. directed edge has (-1) instead of (1)
  signed_pdag = pdag_unsigned_to_signed(cpdag);

  % find all corresponding DAGs
  [n_dags,dag_list] = pdag_to_all_dags( signed_pdag );

