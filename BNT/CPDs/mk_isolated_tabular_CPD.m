function CPD = mk_isolated_tabular_CPD(fam_sz, args)
% function CPD = mk_isolated_tabular_CPD(fam_sz, args)
% function CPD = mk_isolated_tabular_CPD(fam_sz, args)
% Make a single CPD by creating a mini-bnet containing just this one family.
% This is necessary because the CPD constructor requires a bnet.

n = length(fam_sz);
dag = zeros(n,n);
ps = 1:(n-1);
if ~isempty(ps)
  dag(ps,n) = 1;
end
bnet = mk_bnet(dag, fam_sz);
CPD = tabular_CPD(bnet, n, args{:});
