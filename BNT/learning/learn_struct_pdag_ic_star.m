function [pdag, G] = learn_struct_pdag_ic_star(cond_indep, n, k, varargin)
% LEARN_STRUCT_PDAG_IC_STAR Learn a partially oriented DAG (pattern) with latent 
% variables using the IC* algorithm
% P = learn_struct_pdag_ic_star(cond_indep, n, k, ...)
%
% n is the number of nodes.
% k is an optional upper bound on the fan-in (default: n)
% cond_indep is a boolean function that will be called as follows:
% feval(cond_indep, x, y, S, ...)
% where x and y are nodes, and S is a set of nodes (positive integers),
% and ... are any optional parameters passed to this function.
%
% The output P is an adjacency matrix, in which
% P(i,j) = -1 if there is either a latent variable L such that i <-L-> j 
% OR there is a directed edge from i->j.
% P(i,j) = -2 if there is a marked directed i-*>j edge.
% P(i,j) = P(j,i) = 1 if there is and undirected edge i--j
% P(i,j) = P(j,i) = 2 if there is a latent variable L such that i<-L->j.
%
% The IC* algorithm learns a latent structure associated with a set of observed 
% variables. 
% The latent structure revealed is the projection in which every latent variable is
% 1) a root node
% 2) linked to exactly two observed variables.
% Latent variables in the projection are represented using a bidirectional graph, 
% and thus remain implicit.
%
% See Pearl, "Causality: Models, Reasoning, and Inference", 2000, p52 for more details.
% Written by Tamar Kushnir, 2000

sep = cell(n,n);
ord = 0;
done = 0;
G = ones(n,n);
G = setdiag(G,0);
while ~done
  done = 1;
  [X,Y] = find(G); 
  for i=1:length(X)
    x = X(i); y = Y(i);
    nbrs = mysetdiff(myunion(neighbors(G, x), neighbors(G,y)), [x y]);
    if length(nbrs) >= ord & G(x,y) ~= 0
      done = 0;
      SS = subsets(nbrs, ord, ord); % all subsets of size ord
      for si=1:length(SS)
	S = SS{si};
	if feval(cond_indep, x, y, S, varargin{:})  
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
% where y is a common neighbor to x and z
for i=1:length(X)
  x = X(i);
  y = Y(i);
  Z = find(G(y,:));
  Z = mysetdiff(Z, x);
  for z=Z(:)'
    if G(x,z)==0 & ~ismember(y, sep{x,z}) & ~ismember(y, sep{z,x})
      pdag(x,y) = -1; pdag(y,x) = 0;
      pdag(z,y) = -1; pdag(y,z) = 0;
    end
  end
end

% Convert the minimal pattern to a complete one using the following rules:
% Rule 1:
% if a and b are non-adjacent nodes with a common neighbor c,
% if a->c and not b->c then c-*>b (marked arrow).
% Rule 2:
% if a and b are adjacent and there is a directed path (marked links) from a to b
% then a->b (add arrowhead).
%Pearl (2000)

arrowin = [-1 -2 2];
old_pdag = zeros(n);
iter = 0;
while ~isequal(pdag, old_pdag)
  iter = iter + 1;
  old_pdag = pdag;
  % rule 1
  [X, Y] = find(pdag);
  for i=1:length(X)
    x = X(i);
    y = Y(i);
    Z = find(pdag(y,:));
    Z = mysetdiff(Z, x);
    for z=Z(:)'
      if G(x,z)==0 & ismember(pdag(x,y),arrowin) & ~ismember(pdag(z,y),arrowin)
        pdag(y,z) = -2; pdag(z,y) = 0;
      end
    end
  end
  % rule 2
  [X, Y] = find(G); 
  %check all adjacent nodes because if pdag(x,y) = -1 
  %and pdag(y,x) = 0 there could still be an bidirected edge between x & y.
  for i=1:length(X)
    x = X(i);
    y = Y(i);
    if ~ismember(pdag(x,y), arrowin) %x->y doesn't exist yet
      %find marked path from x to y
      add_arrow = marked_path(x,y,pdag);
      if add_arrow 
        if pdag(y,x)==-1 %bidirected edge
          pdag(x,y) = 2; pdag(y,x) = 2;
        else
          pdag(x,y) = -1;pdag(y,x) = 0;
        end
      end
    end
  end
end


%%%%%%%%%%%%%

function t = marked_path(x,y,L)
% MARKED_PATH is a boolean function which returns 1 if a marked path 
% between nodes x and y exists in the partially directed latent structure L.
%
% t = marked_path(x,y,L)
%
% x and y are the starting and ending nodes in the path, respectively.
% L is a latent structure (partially directed graph with possible latent variables).
%
% Rule 2 of IC* algorithm (see Pearl, 2000)

t=0;

%find set of marked links from x
marked = find(L(x,:)==-2);
if ismember(y,marked)
  t=1; %marked path found
else
  for m=marked(:)'
    t = marked_path(m,y,L);
    if t==1
      break; %stop when marked path found
    end
  end
end
