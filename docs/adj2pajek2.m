% ADJ2PAJEK2 Converts an adjacency matrix representation to a Pajek .net read format
% adj2pajek2(adj, filename-stem, 'argname1', argval1, ...)
%
% Set A(i,j)=-1 to get a dotted line
%
% Optional arguments
%
%  nodeNames - cell array, defaults to {'v1','v2,...}
%  shapes - cell array, defaults to {'ellipse','ellipse',...}
%     Choices are 'ellipse', 'box', 'diamond', 'triangle', 'cross', 'empty'
%  partition - vector of integers, defaults to [1 1 ... 1]
%    This will automatically color-code the vertices by their partition
%
% Run pajek (available from http://vlado.fmf.uni-lj.si/pub/networks/pajek/)
% Choose File->Network->Read from the menu
% Then press ctrl-G (Draw->Draw)
% Optional: additionally load the partition file then press ctrl-P (Draw->partition)
%
% Examples
% A=zeros(5,5);A(1,2)=-1;A(2,1)=-1;A(1,[3 4])=1;A(2,5)=1;
% adj2pajek2(A,'foo') % makes foo.net
%
% adj2pajek2(A,'foo','partition',[1 1 2 2 2]) % makes foo.net and foo.clu
% 
% adj2pajek2(A,'foo',...
%            'nodeNames',{'TF1','TF2','G1','G2','G3'},...
%            'shapes',{'box','box','ellipse','ellipse','ellipse'});
%
%
% The file format is documented on p68 of the pajek manual
% and good examples are on p58, p72
%
% Written by Kevin Murphy, 30 May 2007
% Based on adj2pajek by Gergana Bounova
% http://stuff.mit.edu/people/gerganaa/www/matlab/routines.html
% Fixes a small bug (opens files as 'wt' instead of 'w' so it works in windows)
% Also, simplified her code and added some features.

function []=adj2pajek2(adj,filename, varargin)

N = length(adj);
for i=1:N
  nodeNames{i} = strcat('"v',num2str(i),'"');
  shapes{i} = 'ellipse';
end

[nodeNames, shapes, partition] = process_options(varargin, ...
    'nodeNames', nodeNames, 'shapes', shapes, 'partition', []);

if ~isempty(partition)
  fid = fopen(sprintf('%s.clu', filename),'wt','native'); 
  fprintf(fid,'*Vertices  %6i\n',N);
  for i=1:N
    fprintf(fid, '%d\n', partition(i));
  end
  fclose(fid);
end

fid = fopen(sprintf('%s.net', filename),'wt','native'); 

fprintf(fid,'*Vertices  %6i\n',N);
for i=1:N
  fprintf(fid,'%3i %s %s\n', i, nodeNames{i}, shapes{i});
end

%fprintf(fid,'*Edges\n');
fprintf(fid,'*Arcs\n'); % directed
for i=1:N
  for j=1:N
    if adj(i,j) ~= 0
      fprintf(fid,' %4i   %4i   %2i\n',i,j,adj(i,j));
    end
  end
end
fclose(fid)


if 0
adj2pajek2(A,'foo',...
            'nodeNames',{'TF1','TF2','G1','G2','G3'},...
            'shapes',{'box','box','ellipse','ellipse','ellipse'});

N = 100; part = ones(1,N); part(intersect(reg.tfidxTest,1:N))=2;
G = reg.Atest(1:N, 1:N)';
adj2pajek2(G, 'Ecoli100', 'partition', part)
end
