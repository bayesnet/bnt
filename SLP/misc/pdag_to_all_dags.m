function [n_dags,dag_list] = pdag_to_all_dags( pdag )

%
% [n_dags,dag_list] = pdag_to_all_dags( pdag)
%
% generates a cell array of ALL Markov-equivalent DAGs 
% corresponding to a partially directed acyclic graph (PDAG).
%
% Input: PDAG  (PDAG does NOT have to be complete)
%     Format of pdag matrix:
%     Edge with known direction a->b  represented as pdag(a,b)=-1  pdag(b,a)=0
%     Edge with unknown direction a-b represented as pdag(a,b)=1   pdag(b,a)=1
%
% Output: Number of DAGs generated and
%         Cell array of all permissible extensions of PDAG
%
% Sample Use:  
%      % Use output of PC algorithm
%      dag = mk_rnd_dag(4);   % create random DAG
%      % Generate pdag through PC algorithm
%      pdag = learn_struct_pdag_pc('dsep', length(dag), 3, dag)
%      [n_dags,dag_list] = pdag_to_all_dags(pdag);
%      n_dags  
%
% If you want to generate all DAGs that are Markov equivalent to an 
% input DAG (not a pattern), then use function Markov_equivalent_dags(dag)
% instead which calls this function.
%
% =======================================================================
% Algorithm to generate ALL DAGs (pdag_to_all_dags):
%
% 0) Initialize an empty list of DAGs.
%
% 1) Complete current PDAG as far as possible using Rules R1-R4.
%
% 2) Select an unoriented edge X-Y.  
%
%    a) If none left:  
%       Done. Add DAG=abs(PDAG) to list of output DAGs.  Return.
%
%    b) Otherwise:       
%       Select an unoriented edge X-Y.
%       Create PDAG1 with X->Y.
%       Create PDAG2 with Y->X.
%       Recursion: For EACH PDAG (PDAG1/2):  Go to Step 1.
%
%
% This algorithm is a slight modification of the algorithm by Meek (1995)  
% which generates a single DAG extension of a PDAG - here we just add 
% recursion to consider both possible orientations for each considered edge.
%
% For the original algorithm by Meek, see 
%    C. Meek, "Causal inference and causal explanation with background 
%    knowledge", UAI 1995, Section 3.1.1, "Phase III" algorithm.
%
% Thanks to Daniel Eaton for extensive testing and bug reports.
%
% Imme Ebert-Uphoff (ebert@tree.com), 2007
% =======================================================================

   % MAIN
   dag_list={};  % init empty list of DAGs

   % Complete pdag as far as possible using Rules R1-R4 of Meek (1995)
   cpdag = complete_pattern(pdag);
 
  % Start recursion
   dag_list = recurse_unoriented_edge(cpdag,dag_list); 

   % return # of DAGs along with dag_list
   n_dags = length(dag_list);
   if (n_dags == 0)  % no DAGs generated
      fprintf('\nPDAG does not have any permissible extension!\n');
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECURSE_UNORIENTED_EDGE                                        %
%    implements Step 2 of the pdag_to_all_dags algorithm.        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updated_list = recurse_unoriented_edge(cpdag, dag_list)

   % input must be a COMPLETE pdag

   [A,B] = find(cpdag==1);  % find all undirected edges

   updated_list = dag_list;

   if isempty(A)   % if no undirected edges left
      % End of recursion reached.
      % Convert all (-1) values to (1) to yield standard DAG, add DAG to list.
      updated_list{end+1} = abs(cpdag);   

   else
      a = A(1); b = B(1); % choose first unoriented edge 
                          % (any unoriented edge could be used here)

      % choose two different directions for edge and complete BOTH !      
      % PDAG1: contains a -> b      
      %fprintf('PDAG1: %d -> %d\n',a,b);
      pdag1 = cpdag;  
      pdag1(a,b) = -1;  pdag1(b,a) = 0;
      % complete as far as possible using rules R1-R4:
      cpdag1 = complete_pattern(pdag1);  
      % Continue recursion on another unoriented edge
      updated_list = recurse_unoriented_edge(cpdag1,updated_list);

      % PDAG1: contains b -> a
      %fprintf('PDAG2: %d -> %d\n',b,a);
      pdag2 = cpdag;  
      pdag2(a,b) = 0;  pdag2(b,a) = -1;
      % complete as far as possible using rules R1-R4:
      cpdag2 = complete_pattern(pdag2); 
      % Continue recursion on another unoriented edge
      updated_list = recurse_unoriented_edge(cpdag2, updated_list);

   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

