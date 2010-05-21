function cache = score_init_cache(N,L)
% SCORE_INIT_CACHE generate an empty cache for local computation in structure learning
% cache = score_init_cache(number_of_nodes,cache_size)
%
% For 2 nodes with cache of size 5 :
%
% cache =
%   Nw  b        0      0      0 --> Nw=number of writings in cache (+1) and b==1 iff the cache is full
%   0   0        1   -239.12   1 --> 1st familly in the cache (node 1 without parents) calculate with bic
%   0   0        2   -318.98   1
%   1   0        2   -189.23   2 --> 3rd familly in the cache (node 2 with 1 as parent) calculate with bayesian
%   0   1        1   -251.09   1
%   0   0        0      0      0 --> empty entry
%   |   |        |      |      |
%   |   |        |      |      |___> scoring function : 1 for 'bic', 2 for 'bayesian', ...
%   |   |        |      |__________> local score of the familly
%   |   |        |_________________> son node of the familly
%   |   |__________________________> ==1 iff node 2 is parent of son node
%   |______________________________> ==1 iff node 1 is parent of son node
%
%
% V1.1 : 6 may 2003 (O. Francois - francois.olivier.c.h@gmail.com, Ph. Leray - philippe.leray@univ-nantes.fr)
%
%

cache=zeros(L+1,N+3);
cache(1,1)=2;

% using a sparse matrix does not improve performances
