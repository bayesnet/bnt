function c = get_cpts(bnet)
% Get all the cpts in tabular form

cpds = bnet.CPD;
c = cell(size(cpds));
for i = 1:length(c)
  c{i} = CPT(bnet, i);
end
