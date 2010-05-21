function [pdag, G] = learn_struct_pdag_pc_mod(cond_indep, n, k, varargin)
% LEARN_STRUCT_PDAG_PC Learn a partially oriented DAG (pattern) using the PC algorithm
% P = learn_struct_pdag_pc(cond_indep, n, k, ...)
%
% n is the number of nodes.
% k is an optional upper bound on the fan-in (default: n)
% cond_indep is a boolean function that will be called as follows:
%   feval(cond_indep, x, y, S, ...)
% where x and y are nodes, and S is a set of nodes (positive integers),
% and ... are any optional parameters passed to this function.
%
% The output P is an adjacency matrix, in which
% P(i,j) = -1 if there is an i->j edge.
% P(i,j) = P(j,i) = 1 if there is an undirected edge i <-> j
%
% The PC algorithm does structure learning assuming all variables are observed.
% See Spirtes, Glymour and Scheines, "Causation, Prediction and Search", 1993, p117.
% This algorithm may take O(n^k) time if there are n variables and k is the max fan-in,
% but this is quicker than the Verma-Pearl IC algorithm, which is always O(n^n).


sep = cell(n,n);
ord = 0;
done = 0;
G = ones(n,n);
G=setdiag(G,0);
while ~done
 done = 1;
 [X,Y] = find(G);
 for i=1:length(X)
   x = X(i); y = Y(i);
   %nbrs = mysetdiff(myunion(neighbors(G, x), neighbors(G,y)), [x y]);
   nbrs = mysetdiff(neighbors(G, y), x);  % bug fix by Raanan Yehezkel <raanany@ee.bgu.ac.il> 6/27/04
   if length(nbrs) >= ord & G(x,y) ~= 0
     done = 0;
     %SS = subsets(nbrs, ord, ord); % all subsets of size ord
     SS = subsets1(nbrs, ord);
     for si=1:length(SS)
       S = SS{si};
       if feval(cond_indep, x, y, S, varargin{:})
         %if isempty(S)
         %  fprintf('%d indep of %d ', x, y);
         %else
         %  fprintf('%d indep of %d given ', x, y); fprintf('%d ', S);
         %end
         %fprintf('\n');

         % diagnostic
         %[CI, r] = cond_indep_fisher_z(x, y, S, varargin{:});
         %fprintf(': r = %6.4f\n', r);

         G(x,y) = 0;
         G(y,x) = 0;
         sep{x,y} = myunion(sep{x,y}, S);
         sep{y,x} = myunion(sep{y,x}, S);
         break; % no need to check any more subsets
       end
     end
   end
 end
 ord = ord + 1;
end

% Create the minimal pattern,
% i.e., the only directed edges are V structures.
pdag = G;
[X, Y] = find(G);
% We want to generate all unique triples x,y,z
% This code generates x,y,z and z,y,x.
for i=1:length(X)
 x = X(i);
 y = Y(i);
 Z = find(G(y,:));
 Z = mysetdiff(Z, x);
 for z=Z(:)'
   if G(x,z)==0 & ~ismember(y, sep{x,z}) & ~ismember(y, sep{z,x})
     %fprintf('%d -> %d <- %d\n', x, y, z);
     pdag(x,y) = -1; pdag(y,x) = 0;
     pdag(z,y) = -1; pdag(y,z) = 0;
   end
 end
end

% Convert the minimal pattern to a complete one,
% i.e., every directed edge in P is compelled
% (must be directed in all Markov equivalent models),
% and every undirected edge in P is reversible.
% We use the rules of Pearl (2000) p51 (derived in Meek (1995))

old_pdag = zeros(n);
iter = 0;
while ~isequal(pdag, old_pdag)
 iter = iter + 1;
 old_pdag = pdag;
 % rule 1
 [A,B] = find(pdag==-1); % a -> b
 for i=1:length(A)
   a = A(i); b = B(i);
   C = find(pdag(b,:)==1 & G(a,:)==0); % all nodes adj to b but not a
   if ~isempty(C)
     pdag(b,C) = -1; pdag(C,b) = 0;
     %fprintf('rule 1: a=%d->b=%d and b=%d-c=%d implies %d->%d\n', a, b, b, C, b, C);
   end
 end
 % rule 2
 [A,B] = find(pdag==1); % unoriented a-b edge
 for i=1:length(A)
   a = A(i); b = B(i);
   if any( (pdag(a,:)==-1) & (pdag(:,b)==-1)' );
     pdag(a,b) = -1; pdag(b,a) = 0;
     %fprintf('rule 2: %d -> %d\n', a, b);
   end
 end
 % rule 3
 [A,B] = find(pdag==1); % a-b
 for i=1:length(A)
   a = A(i); b = B(i);
   % Bug fix by Imme Ebert-Uphoff (ebert@tree.com), Jan 2007
   % C = find( (G(a,:)==1) & (pdag(:,b)==-1)' );
   C = find( (pdag(a,:)==1) & (pdag(:,b)==-1)' );
   % C contains nodes c s.t. a-c->ba
   G2 = setdiag(G(C, C), 1);
   if any(G2(:)==0) % there are 2 different non adjacent elements of C
     pdag(a,b) = -1; pdag(b,a) = 0;
     %fprintf('rule 3: %d -> %d\n', a, b);
   end
 end
end


%  % Test Rule 3 of PC algorithm
%  
%  % Define PDAG
%  pdag = zeros(4);
%  pdag(2,1)=-1;
%  pdag(3,1)=-1;
%  pdag(2,4)=-1;
%  pdag(3,4)=-1;
%  pdag(1,4)=1;
%  pdag(4,1)=1;
%  
%  fprintf('\nSample input PDAG:\n');
%  pdag
%  
%  fprintf('Sample DAG generated from PDAG:\n');
%  dag = abs(pdag_to_dag(pdag))
%  
%  fprintf('Output from current PC algorithm:\n');
%  pdag_PC = learn_struct_pdag_pc('dsep', 4, 3, dag)
%  
%  % Problem can be fixed by changing Line 120 of learn_struct_pdag_pc.m
%  %    C = find( (G(a,:)==1) & (pdag(:,b)==-1)' );
%  % to
%  %    C = find( (pdag(a,:)==1) & (pdag(:,b)==-1)' );
%  
%  fprintf('Correct version:\n');
%  pdag_PC_mod = learn_struct_pdag_pc_mod('dsep', 4, 3, dag)

