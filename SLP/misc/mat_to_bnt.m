function [res] = mat_to_bnt(mat,misv)
% MAT_TO_BNT Convert a matrix to a cell array
% D = mat_to_bnt(data,misv)
%
% Input : 
%   data(i,m) is the node i in the case m,
%   misv is the way you choose to encode missing data
%               in the original matrix (-9999 by default)
%
% Output :
%  D = cell array containing data(i,m) is the data is OK 
%                or [] is the data is missing
%
%
% V1.2 : 18 feb 2003 (Ph. Leray - philippe.leray@univ-nantes.fr)
%
% >> m=rand(2,4)
%
%   m =
%
%    0.9525    0.4693    0.3907    0.1496
%    0.9274    0.3157    0.1346    0.9383
%
%  >> m(2,2)=-9
%
%   m =
%
%    0.9525    0.4693    0.3907    0.1496
%    0.9274   -9.0000    0.1346    0.9383
%
%  >> mat_to_bnt(m,-9)
%
%   ans = 
%
%    [0.9525]    [0.4693]    [0.3907]    [0.1496]
%    [0.9274]          []    [0.1346]    [0.9383]
%

if nargin <1
    error('Requires at least 1 argument.')
end

if nargin == 1
	misv=-9999;
    end;

taille=size(mat);
long=taille(1);
larg=taille(2);
for i=1:long
  for j=1:larg
    res{i,j}=mat(i,j);
    if(mat(i,j)==misv)
      res{i,j}=[];
    end
  end
end

