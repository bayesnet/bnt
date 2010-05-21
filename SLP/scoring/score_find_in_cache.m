function [bool, score] = score_find_in_cache(cache,j,ps,scoring_fn)
% [bool, score] = score_find_in_cache(cache,j,ps,scoring_fn)
% 
% francois.olivier.c.h@gmail.com

%tic
L=size(cache,1);
N=size(cache,2)-3;

if N<1
  bool=0;
  score=0;
  return
end

parents=zeros(1,N+1);
parents(ps)=1;parents(N+1)=j;

switch scoring_fn
  case 'bic',
    fn=1;
  case 'bayesian',
    fn=2;
  otherwise,
    fn=3;
    %error(['unrecognized scoring fn ' scoring_fn]);     
end

%parent = str2num(num2str(parents,'%1d'));
%[tmp y]=find(cache(:,N+3)==fn);
%if ~isempty(tmp)
%  [tmp2 y]=find(str2num(num2str(cache(tmp,1:N+1),'%1d'))==parent);
%  candidats=tmp(tmp2);
%else
%  candidats=[];
%end

[tmp y]=find(cache(2:L,N+3)==fn);
tmp=tmp+1;
[tmp2 y]=find(cache(tmp,N+1)==j);
candidats=tmp(tmp2);
if ~isempty(candidats)
  for i=1:N      % N=size(cache,2)-3;
    if ~isempty(candidats)
      [tmp2 y]=find(cache(tmp,i)==parents(i));
      candidats=intersect(candidats,tmp(tmp2));
    end
  end
end

%Tpre=toc

if ~isempty(candidats)
  bool=1;
else
  bool=0;
end

if bool
  score=cache(candidats(1),N+2);
else
  score=0;
end