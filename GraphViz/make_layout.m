function [x, y] = layout_dag(adj)
% MAKE_LAYOUT		Creates a layout from an adjacency matrix
%
%  [X, Y] = MAKE_LAYOUT(ADJ)
%
% Inputs :
%	ADJ = adjacency matrix (source, sink)
%
% Outputs :
%	X, Y : Positions of nodes
%
% Usage Example : [X, Y] = make_layout(adj);
%
%
% Note	: Uses some very simple heuristics, so any other
%         algorithm would create a nicer layout 
%
% See also 

% Uses :

% Change History :
% Date		Time		Prog	Note
% 13-Apr-2000	 8:25 PM	ATC	Created under MATLAB 5.3.1.29215a (R11.1)

% ATC = Ali Taylan Cemgil,
% SNN - University of Nijmegen, Department of Medical Physics and Biophysics
% e-mail : cemgil@mbfys.kun.nl 

N = size(adj,1);
tps = toposort(adj);

if ~isempty(tps), % is directed ?
  level = zeros(1,N);
  for i=tps,
    idx = find(adj(:,i));
    if ~isempty(idx),
      l = max(level(idx));
      level(i)=l+1;
    end;
  end;
else
  level = poset(adj,1)'-1;  
end;

y = (level+1)./(max(level)+2);
y = 1-y;
x = zeros(size(y));
for i=0:max(level),
  idx = find(level==i);
  offset = (rem(i,2)-0.5)/10;
  x(idx) = (1:length(idx))./(length(idx)+1)+offset;
end;

%%%%%%%

function [depth] = poset(adj, root)
% POSET		Identify a partial ordering among the nodes of a graph
% 
%  [DEPTH] = POSET(ADJ,ROOT)
% 
% Inputs :
%    ADJ : Adjacency Matrix
%    ROOT : Node to start with
% 
% Outputs :
%    DEPTH : Depth of the Node
% 
% Usage Example : [depth] = poset(adj,12);
% 
% 
% Note     : All Nodes must be connected
% See also 

% Uses :

% Change History :
% Date		Time		Prog	Note
% 17-Jun-1998	12:01 PM	ATC	Created under MATLAB 5.1.0.421

% ATC = Ali Taylan Cemgil,
% SNN - University of Nijmegen, Department of Medical Physics and Biophysics
% e-mail : cemgil@mbfys.kun.nl 

adj = adj+adj';

N = size(adj,1);
depth = zeros(N,1);
depth(root) = 1;
queue = root;

while 1,
  if isempty(queue),
    if all(depth), break; 
    else
      root = find(depth==0); 
      root = root(1);
      depth(root) = 1;
      queue = root;
    end;
  end;
  r = queue(1); queue(1) = [];
  idx = find(adj(r,:));
  idx2 = find(~depth(idx));
  idx = idx(idx2);
  queue = [queue idx];
  depth(idx) = depth(r)+1;
end;

%%%%%%%%%

function [seq] = toposort(adj)
% TOPOSORT		A Topological ordering of nodes in a directed graph
% 
%  [SEQ] = TOPOSORT(ADJ)
% 
% Inputs :
%    ADJ : Adjacency Matrix. 
%	   ADJ(i,j)==1 ==> there exists a directed edge
%	   from i to j
% 
% Outputs :
%    SEQ : A topological ordered sequence of nodes.
%          empty matrix if graph contains cycles.
%
% Usage Example : 
%		N=5;
%		[l,u] = lu(rand(N));
%		adj = ~diag(ones(1,N)) & u>0.5;		
%		seq = toposort(adj);
% 
% 
% Note     :
% See also 

% Uses :

% Change History :
% Date		Time		Prog	Note
% 18-May-1998	 4:44 PM	ATC	Created under MATLAB 5.1.0.421

% ATC = Ali Taylan Cemgil,
% SNN - University of Nijmegen, Department of Medical Physics and Biophysics
% e-mail : cemgil@mbfys.kun.nl 
 
N = size(adj);
indeg = sum(adj,1);
outdeg = sum(adj,2);
seq = [];

for i=1:N,
  % Find nodes with indegree 0
  idx = find(indeg==0);
  % If can't find than graph contains a cycle
  if isempty(idx), 
    seq = [];
    break;
  end;
  % Remove the node with the max number of connections
  [dummy idx2] = max(outdeg(idx));
  indx = idx(idx2);
  seq = [seq, indx];
  indeg(indx)=-1;
  idx = find(adj(indx,:));
  indeg(idx) = indeg(idx)-1;
end;




