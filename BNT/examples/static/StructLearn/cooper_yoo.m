% Do the example in Cooper and Yoo, "Causal discovery from a mixture of experimental and
% observational data", UAI 99, p120

N = 2;
dag = zeros(N);
A = 1; B = 2;
dag(A,B) = 1;
ns = 2*ones(1,N);

bnet0 = mk_bnet(dag, ns);
%bnet0.CPD{A} = tabular_CPD(bnet0, A, 'unif', 1);
bnet0.CPD{A} = tabular_CPD(bnet0, A, 'CPT', 'unif', 'prior_type', 'dirichlet');
bnet0.CPD{B} = tabular_CPD(bnet0, B, 'CPT', 'unif', 'prior_type', 'dirichlet');

samples = [2 2;
	   2 1; 
	   2 2;
	   1 1;
	   1 2;
	   2 2;
	   1 1;
	   2 2;
	   1 2;
	   2 1;
	   1 1];

clamped = [0 0;
	   0 0;
	   0 0;
	   0 0;
	   0 0;
	   1 0;
	   1 0;
	   0 1;
	   0 1;
	   0 1;
	   0 1];

nsamples = size(samples, 1);

% sequential version
LL = 0;
bnet = bnet0;
for l=1:nsamples
  ev = num2cell(samples(l,:)');
  manip = find(clamped(l,:)');
  LL = LL + log_marg_lik_complete(bnet, ev, manip);
  bnet = bayes_update_params(bnet, ev, manip);
end
assert(approxeq(exp(LL), 5.97e-7)) % compare with result from UAI paper


% batch version
cases = num2cell(samples');
LL2 = log_marg_lik_complete(bnet0, cases, clamped');
bnet2 = bayes_update_params(bnet0, cases, clamped');

assert(approxeq(LL, LL2))

for j=1:N
  s1 = struct(bnet.CPD{j}); % violate object privacy
  s2 = struct(bnet2.CPD{j});
  assert(approxeq(s1.CPT, s2.CPT))
end

