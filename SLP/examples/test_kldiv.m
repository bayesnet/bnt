clear all;

b1=mk_asia_bnet;
b2=mk_asia_bnet;
res1=[]; res2=[];

fprintf('\nComputation of several KL_divergences (quick methode 1): ');
tmp=cputime;
for ii=0.1:.1:1,
b2.CPD{1} = tabular_CPD(b2, 1, [ii 1-ii]);
res1=[res1 kl_divergence(b1,b2)];
end
fprintf('%3.2f sec.\n', cputime-tmp);

fprintf('Computation of several KL_divergences (slow methode 2): ');
tmp=cputime;
for ii=0.1:.1:1,
b2.CPD{1} = tabular_CPD(b2, 1, [ii 1-ii]);
res2=[res2 kl_divergence2(b1,b2)];
end
fprintf('%3.2f sec.\n', cputime-tmp);
[res1' res2']