function unprocessed = find_unprocessed(dag)
% unprocessed = find_unprocessed(dag)
%
% francois.olivier.c.h@gmail.com

undirected_edges = dag.*dag';
[unprocessed, tmp] = find(undirected_edges)
unprocessed = unique(unprocessed);
%  N = size(dag,1);
%  unprocessed = [];
%  for i=1:(N-1)
%      for j=(i+1):N
%          if dag(i,j)==1 & dag(j,i)==1
%             unprocessed = [unprocessed,i,j];
%          end
%      end
%  end
%  unprocessed = unique(unprocessed);
