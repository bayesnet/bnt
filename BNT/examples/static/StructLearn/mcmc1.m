% We compare MCMC structure learning with exhaustive enumeration of all dags.

N = 3;
%N = 4;
dag = mk_rnd_dag(N);
ns = 2*ones(1,N);
bnet = mk_bnet(dag, ns);
for i=1:N
  bnet.CPD{i} = tabular_CPD(bnet, i);
end

ncases = 100;
data = zeros(N, ncases);
for m=1:ncases
  data(:,m) = cell2num(sample_bnet(bnet));
end

dags = mk_all_dags(N);
score = score_dags(data, ns, dags);
post  = normalise(exp(score));

[sampled_graphs, accept_ratio] = learn_struct_mcmc(data, ns, 'nsamples', 100, 'burnin', 10);
mcmc_post = mcmc_sample_to_hist(sampled_graphs, dags);

if 0
  subplot(2,1,1)
  bar(post)
  subplot(2,1,2)
  bar(mcmc_post)
  print(gcf, '-djpeg', '/home/cs/murphyk/public_html/Bayes/Figures/mcmc_post.jpg')

  clf
  plot(accept_ratio)
  print(gcf, '-djpeg', '/home/cs/murphyk/public_html/Bayes/Figures/mcmc_accept.jpg')
end
