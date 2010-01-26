function [x, y, h] = draw_dbn(adj, inter, flip_intra, K, labels, node_t, x, y)
% DRAW_LAYOUT_DBN		Draws a layout for a Dynamical Belief Network
%
%  [<X1, Y1, X2, Y2>] = DRAW_LAYOUT_DBN(INTRA, INTER, <FLIP_FLAG, K, LABELS, ISBOX, X1, Y1>)
%
% Inputs :
%	INTRA, INTER : Adjacency matrices
%       FLIP_FLAG : Transposes the DAG layout obtained from INTRA connections
%                   If X1, Y1 are specified, FLIP_FLAG has no effect.
%       K         : Unfold K times <Default = 2>
%       LABELS - if -1, we use 1:N*K
%       Rest : See DRAW_LAYOUT      
% 
% Outputs :
%	 Xi, Yi : Coordinates of nodes (for i'th timeslice) on the unit square
%        H      : Object Handles
%
% Usage Example : draw_layout_dbn(intra, inter, 1);
%                 draw_layout_dbn(intra, inter);
%
% Note	:
% See also DRAW_GRAPH

% Uses : DRAW_GRAPH

% Change History :
% Date		Time		Prog	Note
% 17-Apr-2000	 1:02 PM	ATC	Created under MATLAB 5.3.1.29215a (R11.1)

% ATC = Ali Taylan Cemgil,
% SNN - University of Nijmegen, Department of Medical Physics and Biophysics
% e-mail : cemgil@mbfys.kun.nl 

N = size(adj,1);
if nargin<3, 
  flip_intra = 0;
end;

if nargin<4,
  K = 2;
end;

if K<2 | K>7, error('2<=K<=7 must hold..'); end;


if nargin<5 
%  labels = cellstr(char(zeros(N,1)+double('+')));
%  labels = cellstr(int2str((1:N)'));
    labels = cellstr(char((0:N-1)'+double('a')));
end;

if nargin<6,
  node_t = zeros(N,1);
%  node_t = rand(N,1) > 0.5;
end;

if nargin<7,
  [x1 y1] = make_layout(adj);
  if flip_intra, tmp = x1; x1 = y1; y1 = tmp; end;
end;

mid = round(K/2);


xi = x1(:)-1;
x = [];
y = repmat(y1(:), [K 1]);
node_t2 = repmat(node_t(:), [K 1]);

if isa(labels,'double') & labels==-1 % KPM
  lb = num2strcell(1:N*K);
else
  lb = {};
  for i=1:K,
    labels1 = labels(:);
    if i==mid,     str = ''; else str = sprintf('%+d',i-mid); end;
    for i=1:N,
      labels1{i} = [labels1{i} '(t' str ')'];  
    end;
    lb = [lb; labels1(:)];
  end;
end

dag = zeros(N*K);

for i=1:K,
  xi = xi+1;
  x = [x; xi];

  idx = ((i-1)*N+1):i*N;
  dag(idx,idx) = adj;
  if i<K,
    dag(idx,idx+N) = inter;
  end;
end;

[x, y, h] = draw_graph(dag, lb, node_t2, x/K, y);


