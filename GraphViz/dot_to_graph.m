function [Adj, labels, x, y] = dot_to_graph(filename)
% [Adj, labels, x, y] = dot_to_graph(filename)
% Extract a matrix representation, node labels, and node position coordinates
% from a file in GraphViz format http://www.research.att.com/sw/tools/graphviz
%
% INPUTS:
%    'filename' - the file in DOT format containing the graph layout.
% OUTPUT:
%  'Adj'    - an adjacency matrix representation of the graph in 'filename';
% 'labels'  - a character array with the names of the nodes of the graph;
%    'x'    - a row vector with the x-coordinates of the nodes in 'filename';
%    'y'    - a row vector with the y-coordinates of the nodes in 'filename'.
%
% WARNINGS: not guaranted to parse ANY GraphViz file. Debugged on undirected 
%       sample graphs from GraphViz(Heawood, Petersen, ER, ngk10_4, process). 
%       Complaines about RecursionLimit set only to 500 on huge graphs.
%       Ignores singletons (disjoint nodes).         
% Sample DOT code "ABC.dot", read by [Adj, labels, x, y] = dot_to_graph('ABC.dot')
% digraph G {
%       A [pos="28,31"];
%       B [pos="74,87"];
%       A -- B [pos="e,61,71 41,47 46,53 50,58 55,64"];
% }
%                                                     last modified: Jan 2004
% by Alexi Savov:  asavov @wustl.edu  |  http://artsci.wustl.edu/~azsavov
%    Leon Peshkin: pesha @ai.mit.edu  |  http://www.ai.mit.edu/~pesha 
%    Tom Minka

if ~exist(filename)                % Checks whether the specified file exists.
  error('* * * File does not exist or could not be found. * * *');
end;

lines = textread(filename,'%s','delimiter','\n','commentstyle','c');  % Read file into cell array
dot_lines = strvcat(lines);                                % of lines, ignoring C-style comments

if findstr(dot_lines(1,:), 'graph ') == []           % Is this a DOT file ?
  error('* * * File does not appear to be in valid DOT format. * * *');
end;

Nlns = size(dot_lines,1);             % The number of lines;
nodes = {};
unread = 1:Nlns;             % 'unread' list of lines which has not been examined yet
edge_id = 1;
Adj = [];
for line_ndx = 1:Nlns   % This section sets the adjacency matrix A(Lnode,Rnode) = edge_id.
  line = dot_lines(line_ndx,:);
  Ddash_pos = strfind(line, ' -- ') + 1;    % double dash positions
  arrow_pos = strfind(line, ' -> ') + 1;    % arrow  dash positions
  tokens = strread(line,'%s','delimiter',' "');
  left_bound = 1;
  for dash_pos = [Ddash_pos arrow_pos];  % if empty - not a POS line
    Lnode = sscanf(line(left_bound:dash_pos -2), '%s');
    Rnode = sscanf(line(dash_pos +3 : length(line)-1),'%s',1);
    Lndx = strmatch(Lnode, nodes, 'exact');
    Rndx = strmatch(Rnode, nodes, 'exact');
    if isempty(Lndx)         % extend our list of nodes 
      nodes{end+1} = Lnode;
      Lndx = length(nodes);
    end
    if isempty(Rndx)
      nodes{end+1} = Rnode;
      Rndx = length(nodes);
    end
    Adj(Lndx, Rndx) = edge_id;
    if  ismember(dash_pos, Ddash_pos)   % The edge is undirected, A(Rndx,LndxL) is also set to 1;
      Adj(Rndx, Lndx) = edge_id;
    end
    edge_id = edge_id + 1; 
    left_bound = dash_pos + 3;
    unread = setdiff(unread, line_ndx); 
  end
end
Nvrt = length(nodes);    % number of vertices we found  [Do we ever have singleton vertices ???]
% nodes = strvcat(nodes); % convert to the searchable array
x = zeros(1, Nvrt); 
y = zeros(1, Nvrt);
labels = nodes;
% Find node's position coordinates if they are contained in 'filename'.
for line_ndx = unread        % Look for node's coordinates among the 'unread' lines.
  line = dot_lines(line_ndx,:);
  bra_pos = strfind(line, '[');       % has to have "[" if it has the label
  lst_node = 0;
  for node = 1:Nvrt     % look through the list of nodes 
    %  THE NEXT STATEMENT we assume no node is substring of any other node
    lbl_pos = strfind(line, nodes{node});
    if (~isempty(lbl_pos) & ~isempty(bra_pos) & (x(node) == 0))  % make sure we have not seen it 
      if (lbl_pos(1) < bra_pos(1))  % label has to be to the left of bracket
	lst_node = node;
      end
    end
  end
  if lst_node
    pos_pos = strfind(line, 'pos');     % position of the "pos"
    if ~isempty(pos_pos)   % this line contains SOME position  
      [node_pos] = sscanf(line(pos_pos:end), ' pos  = "%d,%d"')';
      x(lst_node) = node_pos(1);
      y(lst_node) = node_pos(2);
    end
    % minka
    label_pos = strfind(line, 'label'); % position of the "label"
    if ~isempty(label_pos)
      label_end = strfind(line(label_pos:end),',');
      labels{lst_node} = unquote(line(label_pos+(6:label_end(1)-2)));
    end
  end
end

if (isempty(find(x)) & (nargout > 2))   % If coordinates were requested, but not found in 'filename'.
  warning('File does not contain node coordinates.');
end;
if ~(size(Adj,1)==size(Adj,2))           % Make sure Adj is a square matrix. ? 
  Adj = eye(max(size(Adj)),size(Adj,1))*Adj*eye(size(Adj,2),max(size(Adj)));
end;
x = .9*(x-min(x))/range(x)+.05;  % normalise and push off margins 
y = .9*(y-min(y))/range(y)+.05; 



function s = unquote(s)

s = strrep(s,'"','');
