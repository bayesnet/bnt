function cache = gener_empty_cache(N,L)
% cache = gener_empty_cache(number_of_nodes,lenght_of_cache)
%
% exemple for 2 nodes with cache of size 5 :
%
% cache =
%   5   b        0      0      0 --> 1st empty place and b==1 iff the cache is full
%   0   0        1   -239.12   1 --> 1st familly in the cache (node 1 without parents) calculate with bic
%   0   0        2   -318.98   1
%   1   0        2   -189.23   2 --> 3rd familly in the cache (node 2 with 1 as parent) calculate with bayesian
%   0   1        1   -251.09   1
%   0   0        0      0      0 --> empty entry
%   |   |        |      |      |
%   |   |        |      |      |___> 1 for 'bic' or 2 for 'bayesian'
%   |   |        |      |__________> score of the familly
%   |   |        |_________________> son node of the familly
%   |   |__________________________> ==1 iff node 2 is parent of son node
%   |______________________________> ==1 iff node 1 is parent of son node
%
%
%
%  designed by francois.olivier.c.h@gmail.com
%

cache=zeros(L+1,N+3);
cache(1,1)=2;