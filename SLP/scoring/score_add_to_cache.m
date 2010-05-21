function [cache, place] = score_add_to_cache(cache,j,ps,score,scoring_fn)
% [cache place] = score_add_to_cache(cache,j,ps,score,scoring_fn)
% 
% j is the son node,
% ps is the list of parents of j, for example [12 5 7],
% score is the score to add for this familly.
% scoring_fn is 'bic' or 'bayesian'.
%
% place = where the entry was add.
%
% example for 2 nodes with cache of size 5 :
%
% cache =
%   5   b        0      0      0  --> number of writing in cache +1 and b==1 iff the cache is full
%   0   0        1   -239.12   1  --> 1st familly in the cache (node 1 without parents) calculate with bic
%   0   0        2   -318.98   1
%   1   0        2   -189.23   2  --> 3rd familly in the cache (node 2 with 1 as parent) calculate with bayesian
%   0   1        1   -251.09   1
% .ps2bool.      j    score  1or2 --> new entry
%   |   |        |      |      |
%   |   |        |      |      |___> 1 for 'bic' or 2 for 'bayesian'
%   |   |        |      |__________> score of the familly
%   |   |        |_________________> son node of the familly
%   |   |__________________________> ==1 iff node 2 is parent of son node
%   |______________________________> ==1 iff node 1 is parent of son node
%
% If the cache is FULL then the new place is RanDoMly choose.
%
% francois.olivier.c.h@gmail.com

N=size(cache,2)-3;
place=0;

if ~isempty(find(ps==j))
  disp('This is a cyclic entry, nothing was done.');
elseif j>N | j<0
  disp('This entry is not valid, nothing was done.');
else

  switch scoring_fn
    case 'bic',
      fn=1;
    case 'bayesian',
      fn=2;
    otherwise,
      fn=3;
      %error(['unrecognized scoring fn ' scoring_fn]);
  end
  L=size(cache,1);

  if cache(1,2)==0
    place=cache(1,1);
  else
    place=ceil(rand(1)*(L-1))+1;
  end

  cache(place,:)=0
  cache(place,ps)=1;
  cache(place,N+1)=j;
  cache(place,N+2)=score;
  cache(place,N+3)=fn;
  cache(1,1)=place+1;
  if place==L | cache(1,2)~=0
    cache(1,2)=1
  end

end
