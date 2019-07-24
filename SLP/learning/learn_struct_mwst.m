function [T, score_mat] = learn_struct_mwst(data, discrete, node_sizes, node_type, scoring_fn, root)
% LEARN_STRUCT_MWST Learn an oriented tree using the MSWT algorithm
% T = learn_struct_mwst(data, discrete, node_sizes, node_type, scoring_fn, root)
%
% Input : 
%   data(i,m) is the node i in the case m,
%   discrete = [ 1 if discret-node 0 if not ],
%   node_sizes = 1 if gaussian node,
%   node_type = {'tabular','gaussian',...},
%   score = 'bic' (for complete data and any node types) or 'mutual_info' (tabular nodes),
%   root is the futur root-node of the tree T.
%
% Output :
%	T = adjacency matrix of the tree
%
% V1.2 : 17 feb 2003 (O. Francois - francois.olivier.c.h@gmail.com, Ph. Leray - philippe.leray@univ-nantes.fr)
%
%
% See Chow&Liu 1968 for the original algorithm using Mutual Information scoring.
% Or Heckerman 1994.

if nargin <4
    error('Requires at least 4 arguments.')
end

if nargin == 4
    scoring_fn='bic'; root=1;
end;

if nargin == 5
    root=1;
end;


N=size(data,1);
score_mat=zeros(N,N);

switch scoring_fn
case 'bic',
    for i=1:(N-1)
            score2 = score_family(i, [], node_type{i}, scoring_fn, node_sizes, discrete, data,[]);
        for j=(i+1):N
            score1 = score_family(i, [j], node_type{i}, scoring_fn, node_sizes, discrete, data,[]);
            score = score2-score1;
            score_mat(i,j)=score;
            score_mat(j,i)=score;
        end
    end
case 'mutual_info',
    for i=1:(N-1)
        for j=(i+1):N
            score_mat(i,j)= -mutual_info_score(i,node_sizes(i),j,node_sizes(j),data);
            score_mat(j,i)=score_mat(i,j);
        end
    end
otherwise,
    error(['unrecognized scoring fn ' scoring_fn]);
end

G = minimum_spanning_tree(score_mat);
T = mk_rooted_tree(G, root);
T=full(T);