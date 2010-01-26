function draw_dot(adj);
%
%  draw_dot(name) 
%  
% Sample code illustrating use of dot_to_graph.m function
% Leon Peshkin  
if ispc, shell = 'dos'; else, shell = 'unix'; end  %  Which OS ?

cmdline = strcat(shell,'(''neato -V'')');
status = eval(cmdline);
[status, result] = dos('neato -V');  % request version to check NEATO
if status == 1,  fprintf('Complaining \n'); exit, end

tmpDOTfile = '_GtDout.dot';            % to be platform independant no use of directories
tmpLAYOUT  = '_LAYout.dot'; 
directed = 0;                          % assume UN-directed graph
graph_to_dot(adj > 0, 'directed', directed, 'filename', tmpDOTfile);  % save in file

cmdline = strcat([shell '(''neato -Tdot ' tmpDOTfile ' -o ' tmpLAYOUT ''')']); % preserve trailing spaces 
status = eval(cmdline);         %  get NEATO todo layout

[adj, labels, x, y] = dot_to_graph(tmpLAYOUT);  %  load layout 
delete(tmpLAYOUT); delete(tmpDOTfile);     % clean up temporary files

figure(1); clf; axis square      %  now plot 
[x, y, h] = draw_graph(adj>0, labels, zeros(size(x,2),1), x, y);