function res = bnt_to_mat(data,misv)
% res = bnt_to_mat(data,misv)
%

if nargin <1, error('Requires at least 1 argument.'), end
if nargin == 1,	misv=-9999; end
taille=size(data);
long=taille(1);
larg=taille(2);
for i=1:long
  for j=1:larg
    if ~isempty(data{i,j})
      res(i,j)=data{i,j};
    else
      res(i,j)=misv;
    end
  end
end
