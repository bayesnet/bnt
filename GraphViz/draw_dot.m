function draw_dot(adj,varargin);
%DRAW_DOT   Draw a graph.
% DRAW_DOT(ADJ) plots the graph ADJ in the current figure window, using 
% 'neato' to optimize the layout.
%
% Optional arguments can be passed as name/value pairs: [default]
%
% 'isbox'     - a vector specifying which nodes should be boxed [0]
% 'rotate'    - rotate the graph so that nodes are vertically aligned [1]
% 'tolerance' - alignment tolerance for 'rotate' [0.001]
% 'start'     - a random seed (to select different solutions)
% 'options'   - a string of command-line options for 'neato' ['']
% All of the optional arguments to graph_to_dot are also supported, such as
% 'node_label'.
%
% See also GRAPH_TO_DOT.
%
% Example:
% size=15; Adj = rand(size) > .8;
% Adj2 = triu(Adj,1)+ triu(Adj,1)' + diag(zeros(size,1));
% draw_dot(Adj2)

% Original: Leon Peshkin  
% Modified by Tom Minka

% minka
N = size(adj,1);
unique_labels = cellstr(num2str((1:N)','%-1d'));
labels = unique_labels;
isbox = zeros(N,1);
rotate_flag = 1;
tolerance = 0.001;
options = '';
for i = 1:2:length(varargin)
  switch varargin{i}
    case 'node_label', labels = varargin{i+1}; 
      % replace with unique labels
      varargin{i+1} = unique_labels;
    case 'isbox', isbox = varargin{i+1};
    case 'rotate', rotate_flag = varargin{i+1};
    case 'tolerance', tolerance = varargin{i+1};
    case 'start', start = varargin{i+1}; 
      options = [options ' -Gstart=' num2str(start)];
    case 'options', options = [options ' ' varargin{i+1}];
  end
end

if ispc, shell = 'dos'; else, shell = 'unix'; end  %  Which OS ?

cmdline = strcat(shell,'(''neato -V'')');
status = eval(cmdline);
%[status, result] = dos('neato -V');  % request version to check NEATO
if status == 1,  fprintf('Complaining \n'); exit, end

tmpDOTfile = '_GtDout.dot';            % to be platform independant no use of directories
tmpLAYOUT  = '_LAYout.dot'; 
graph_to_dot(adj > 0, 'filename', tmpDOTfile, 'node_label', unique_labels, varargin{:});  % save in file

cmdline = strcat([shell '(''neato -Tdot ' tmpDOTfile options ' -o ' tmpLAYOUT ''')']); % preserve trailing spaces 
status = eval(cmdline);         %  get NEATO todo layout

[adj, permuted_labels, x, y] = dot_to_graph(tmpLAYOUT);  %  load layout 
delete(tmpLAYOUT); delete(tmpDOTfile);     % clean up temporary files

% permute the original arguments to match permuted_labels.
order = [];
for i = 1:length(permuted_labels)
  j = strmatch(permuted_labels{i},unique_labels,'exact');
  order(i) = j(1);
end
labels = labels(order);
isbox = isbox(order);
if rotate_flag
  [x,y] = best_rotation(x,y,tolerance);
end

figure(1); clf; axis square      %  now plot 
[x, y, h] = draw_graph(adj>0, labels, isbox, x, y, varargin{:});


function [x,y] = best_rotation(x,y,h)
% Rotate the points to maximize the horizontal and vertical alignment.
% Written by Tom Minka.

xm = mean(x);
ym = mean(y);
xr = max(x)-min(x);
yr = max(y)-min(y);
x = (x-xm)/xr;
y = (y-ym)/yr;

xy = [x(:) y(:)];
if 1
  angle = fminbnd(@rotation_cost,-pi/4,pi/4,[],xy,h);
else
  angles = linspace(-pi/4,pi/4,40);
  e = [];
  for i = 1:length(angles)
    e(i) = rotation_cost(angles(i),xy,h);
  end
  %figure(2)
  %plot(angles*180/pi,e)
  angle = angles(argmin(e));
end
%angle*180/pi
c = cos(angle); s = sin(angle);
xy = xy*[c s; -s c];

x = xy(:,1)*xr+xm;
y = xy(:,2)*yr+ym;


function e = rotation_cost(angle,xy,h)
% xy is 2-column matrix.
% e is small if many x's and y's are aligned.

c = cos(angle); s = sin(angle);
xy = xy*[c s; -s c];
dx = sqdist(xy(:,1)',xy(:,1)');
dy = sqdist(xy(:,2)',xy(:,2)');
dx = setdiag(dx,Inf);
dy = setdiag(dy,Inf);
e = sum(exp(-dx(:)/h))+sum(exp(-dy(:)/h));
e = -e;
