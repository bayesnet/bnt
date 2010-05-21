function dag = learn_struct_tan(data, class_node, root, node_sizes, scoring_fn)
% LEARN_STRUCT_TAN Learn the structure of the tree augmented naive bayesian network 
% (with discrete nodes)
% dag = learn_struct_tan(app, class, root, node_sizes)
%
% Input :
% 	data(i,m) is the value of node i in case m
% 	class_node is the class node
% 	root is the root node of the tree part of the dag (must be different from the class node)
%   	node_sizes = 1 if gaussian node,
%   	scoring_fn = 'bic' (default value) or 'mutual_info'
%
% Output :
%	dag = adjacency matrix of the dag
%
% V1.1 : 21 may 2003, (O. Francois - francois.olivier.c.h@gmail.com, Ph. Leray - philippe.leray@univ-nantes.fr)
% V1.2 : may 2005 bug correction about node types (Navid Serrano <Navid.Serrano@jpl.nasa.gov>)


if nargin <4
    error('Requires at least 4 arguments.')
end

if nargin == 4
    scoring_fn='bic';
end;

if class_node==root
    error(' The root node can''t be the class node.');
end

%  if root>class_node
%      root=root-1;
%  end

N=size(data,1);
node_types=cell(N-1,1);
notclass=setdiff(1:N,class_node);
for i=1:N
    if node_sizes(i)==1
        node_types{i}='gaussian';
    else
        node_types{i}='tabular';
    end
end

dag=zeros(N);
T = learn_struct_mwst4tan(data, ones(1,N), node_sizes, node_types, scoring_fn, root, class_node);
dag=T;
dag(class_node,notclass)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T, score_mat] = learn_struct_mwst4tan(data, discrete, node_sizes, node_type, scoring_fn, root, class)

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
score_mat(class,:)=Inf;
score_mat(:,class)=Inf;

switch scoring_fn
case 'bic',
    for i=mysetdiff(1:(N-1), class)
            score2 = score_family(i, [class], node_type{i}, scoring_fn, node_sizes, discrete, data,[]);
        for j=mysetdiff((i+1):N, class)
            score1 = score_family(i, [j,class], node_type{i}, scoring_fn, node_sizes, discrete, data,[]);
            score = score2-score1;
            score_mat(i,j)=score;
            score_mat(j,i)=score;
        end
    end
case 'mutual_info',   % tabular nodes only
    for i=mysetdiff(1:(N-1), class)
        for j=mysetdiff((i+1):N, class)
            score_mat(i,j)= -cond_mutual_info_score(i,node_sizes(i),j,node_sizes(j),class,node_sizes(class),data);
            score_mat(j,i)=score_mat(i,j);
        end
    end
otherwise,
    error(['unrecognized scoring fn ' scoring_fn]);
end

variab = mysetdiff(1:N,class);
%score_mat
G = minimum_spanning_tree(score_mat(variab,variab));
if root>class, root=root-1;end
T = mk_rooted_tree(G, root);
T1 = full(T);
T=zeros(N);
T(variab,variab)=T1;



