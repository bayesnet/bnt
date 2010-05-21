function [n,xechan] = histc_ic(x,edges)
%HISTC_IC  Histogram count 
%
%   [N XECHAN] = HISTC_IC(X,EDGES), for vector X, counts the number of values
%   in X that fall between the elements in the EDGES vector 
%   (EDGES = cell array returned by HIST_IC function)
%
%   N is cell array containing these counts.
%    (or a vector if X is a column vector)
%
%   XECHAN = discretized version of X
%
%   Example :
%	X=randn(100,5);
%	Xapp=X(1:50,:);
%	Xtest=X(51:100,:);
%
%       % bins are computed with Xapp data
%	[n1,bornes,nbbornes,xechan]=hist_ic(Xapp);
%       % histogram is computed on Xtest data, with previously obtained bins
%	[n2,xtechan]=histc_ic(Xtest,bornes);
%
%   05-06-2001 by Ph. Leray - philippe.leray@univ-nantes.fr
%

[nb_l,nb_c]=size(x);
% Outputs declaration
xechan=zeros(nb_l,nb_c);
if nb_c==1, edges={edges}; end
% edges=mat2cell(edges); modified by francois.olivier.c.h@gmail.com
for j=1:nb_c,
  [n{j} xechan(:,j)]=histc(x(:,j),[-inf edges{j} inf]);
  n{j}=n{j}(1:end-1);
end
if nb_c==1, n=n{1}; end
