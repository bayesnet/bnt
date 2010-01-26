function graph_to_dot(G, varargin)
% DAG_TO_DOT Make a file representing the directed graph in dotty format.
% dag_to_dot(G, ...)
%
% Optional arguments should be passed as name/value pairs [default]
%
% 'filename' - if omitted, we write to 'tmp.dot', convert this to 'tmp.ps',
%              and then call ghostview automatically 
% 'arc_label' - arc_label{i,j} is a string attached to the i->j arc. [""]
% 'node_label' - node_label{i} is a string attached to node i. ["i"]
% 'width'      - width in inches [10]
% 'height'     - height in inches [10]
% 'leftright'  - 1 means layout left-to-right, 0 means top-to-bottom [0]
% 'directed'  - 1 means use directed arcs, 0 means undirected [1]
%
% For details on dotty, See http://www.research.att.com/sw/tools/graphviz
%
% Example:
% G = rand(5,5);
% names = cell(5,5);
% names{1,2} = 'arc 1-2';
% graph_to_dot(G, 'arc_label', names)
% or graph_to_dot(G, 'arc_label', 'numbers') % prints value of G(i,j) on i->j arc 

% Kevin Murphy, 1998

% set default args
filename = [];
node_label = [];
arc_label = [];
width = 10;
height = 10;
leftright = 0;
directed = 1;
% get optional args
args = varargin;
for i=1:2:length(args)
  switch args{i}
   case 'filename', filename = args{i+1};
   case 'node_label', node_label = args{i+1};
   case 'arc_label', arc_label = args{i+1};
   case 'width', width = args{i+1};
   case 'height', height = args{i+1};
   case 'leftright', leftright = args{i+1};
   case 'directed', directed = args{i+1};
  end
end

if isstr(arc_label) & strcmp(arc_label, 'numbers')
  N = length(G);
  arc_label = cell(N,N);
  for i=1:N
    for j=1:N
      arc_label{i,j} = sprintf('%4.2f', G(i,j));
    end
  end
end

if isempty(filename)
  make_file(G, 'tmp.dot', node_label, arc_label, width, height, leftright, directed);
  if isunix
    !dot -Tps tmp.dot -o tmp.ps

    !gs tmp.ps &
  else
    dos('dot -Tps tmp.dot -o tmp.ps');
    dos('gsview32 tmp.ps &');
  end
else
  
  
  make_file(G, filename, node_label, arc_label, width, height, leftright, directed);
end


%%%%%%

function make_file(G, filename, node_label, arc_label, width, height, leftright, directed)

n = length(G);
fid = fopen(filename, 'w');
if directed
  fprintf(fid, 'digraph G {\n');
else
  fprintf(fid, 'graph G {\n');
end
fprintf(fid, 'center = 1;\n');
fprintf(fid, 'size=\"%d,%d\";\n', width, height);
if leftright
  fprintf(fid, 'rankdir=LR;\n');
end
for i=1:n
  if isempty(node_label)
    fprintf(fid, '%d;\n', i);
  else
    fprintf(fid, '%d [ label = "%s" ];\n', i, node_label{i});
  end
end
if directed
  for i=1:n
    cs = children(G,i);
    for j=1:length(cs)
      c = cs(j);
      if isempty(arc_label)
	fprintf(fid, '%d -> %d;\n', i, c);
      else
	fprintf(fid, '%d -> %d [label="%s"];\n', i, c, arc_label{i,c});
      end
    end
  end
else
  for i=1:n
    ns = intersect(neighbors(G,i), i+1:n); % remove duplicate arcs
    for j=1:length(ns)
      c = ns(j);
      if isempty(arc_label)
	fprintf(fid, '%d -- %d [dir=none];\n', i, c);
      else
	fprintf(fid, '%d -- %d [label="%s",dir=none];\n', i, c, arc_label{i,c});
      end
    end
  end
end
fprintf(fid, '\n}');
fclose(fid);



%%%%%%%%%%%%%%%

function cs = children(adj_mat, i, t)
% CHILDREN Return the indices of a node's children in sorted order
% c = children(adj_mat, i, t)
%
% t is an optional argument: if present, dag is assumed to be a 2-slice DBN

if nargin < 3 
  cs = find(adj_mat(i,:));
else
  if t==1
    cs = find(adj_mat(i,:));
  else
    ss = length(adj_mat)/2;
    j = i+ss;
    cs = find(adj_mat(j,:)) + (t-2)*ss;
  end
end

%%%%%%%%%%%%

function ps = parents(adj_mat, i)
% PARENTS Return the list of parents of node i
% ps = parents(adj_mat, i)

ps = find(adj_mat(:,i))';

%%%%%%%%%%%%%

function ns = neighbors(adj_mat, i)
% NEIGHBORS Find the parents and children of a node in a graph.
% ns = neighbors(adj_mat, i)

ns = union(children(adj_mat, i), parents(adj_mat, i));



