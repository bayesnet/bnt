function graph_to_dot(adj, varargin)
%GRAPH_TO_DOT  Makes a GraphViz (AT&T) file representing an adjacency matrix
% graph_to_dot(adj, ...) writes to the specified filename.
%
% Optional arguments can be passed as name/value pairs: [default]
%
% 'filename' - if omitted, writes to 'tmp.dot'
% 'arc_label' - arc_label{i,j} is a string attached to the i-j arc [""]
% 'node_label' - node_label{i} is a string attached to the node i ["i"]
% 'width'     - width in inches [10]
% 'height'    - height in inches [10]
% 'leftright' - 1 means layout left-to-right, 0 means top-to-bottom [0]
% 'directed'  - 1 means use directed arcs, 0 means undirected [1]
%
% For details on graphviz, See http://www.research.att.com/sw/tools/graphviz
%
% See also dot_to_graph and draw_dot.

% First version written by Kevin Murphy 2002.
% Modified by Leon Peshkin, Jan 2004.
% Bugfix by Tom Minka, Mar 2004.
                   
node_label = [];   arc_label = [];   % set default args
width = 10;        height = 10;
leftright = 0;     directed = 1;     filename = 'tmp.dot';

for i = 1:2:nargin-1                    % get optional args
  switch varargin{i}
    case 'filename', filename = varargin{i+1};
    case 'node_label', node_label = varargin{i+1};
    case 'arc_label', arc_label = varargin{i+1};
    case 'width', width = varargin{i+1};
    case 'height', height = varargin{i+1};
    case 'leftright', leftright = varargin{i+1};
    case 'directed', directed = varargin{i+1};
  end
end
% minka
if ~directed
  adj = triu(adj | adj');
end

fid = fopen(filename, 'w');
if directed
  fprintf(fid, 'digraph G {\n');
  arctxt = '->'; 
  if isempty(arc_label)
    labeltxt = '';
  else
    labeltxt = '[label="%s"]';
  end
else
  fprintf(fid, 'graph G {\n');
  arctxt = '--'; 
  if isempty(arc_label)
    labeltxt = '[dir=none]';
  else
    labeltext = '[label="%s",dir=none]';
  end
end
edgeformat = strcat(['%d ',arctxt,' %d ',labeltxt,';\n']);
fprintf(fid, 'center = 1;\n');
fprintf(fid, 'size=\"%d,%d\";\n', width, height);
if leftright
  fprintf(fid, 'rankdir=LR;\n');
end
Nnds = length(adj);
for node = 1:Nnds               %  process nodes 
  if isempty(node_label)
    fprintf(fid, '%d;\n', node);
  else
    fprintf(fid, '%d [ label = "%s" ];\n', node, node_label{node});
  end
end
for node1 = 1:Nnds   % process edges
  arcs = find(adj(node1,:));         % children(adj, node);
  for node2 = arcs
    if  ~isempty(arc_label)
      fprintf(fid, edgeformat,node1,node2,arc_label{node1,node2});
    else
      fprintf(fid, edgeformat, node1, node2);    
    end    
  end
end
fprintf(fid, '}');
fclose(fid); 
