function signed_pdag = pdag_unsigned_to_signed(pdag)

% Input: pdag with '1' in one place for every directed edge and 
%                  '1' in two places for every undirected edge
% 
% Output: pdag with '-1' in one place for every directed edge
%                   '1' in two places for every undirected edge
%
% This function is used by Markov_equivalent_dags(dag) to convert 
% output of SLP function dag_to_cpdag(dag) 
% to format required as input for pdag_to_all_dags(pdag).
%
% I'm sure there's a prettier way to code this!
%
% Imme Ebert-Uphoff (ebert@tree.com), 2007
%

  undirected = ( (pdag+pdag')/2 == 1); % extract undirected eges
  directed   = pdag - undirected;       % extract directed edges
  signed_pdag = undirected - directed; % 1s for undirected, (-1)s for directed


