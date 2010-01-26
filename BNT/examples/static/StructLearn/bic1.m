% compare BIC and Bayesian score 

N = 4;
dag = zeros(N,N);
%C = 1; S = 2; R = 3; W = 4; % topological order
C = 4; S = 2; R = 3; W = 1; % arbitrary order
dag(C,[R S]) = 1;
dag(R,W) = 1;
dag(S,W)=1;


false = 1; true = 2;
ns = 2*ones(1,N); % binary nodes
bnet = mk_bnet(dag, ns);
bnet.CPD{C} = tabular_CPD(bnet, C, 'CPT', [0.5 0.5]);
bnet.CPD{R} = tabular_CPD(bnet, R, 'CPT', [0.8 0.2 0.2 0.8]);
bnet.CPD{S} = tabular_CPD(bnet, S, 'CPT', [0.5 0.9 0.5 0.1]);
bnet.CPD{W} = tabular_CPD(bnet, W, 'CPT', [1 0.1 0.1 0.01 0 0.9 0.9 0.99]);


seed = 0;
rand('state', seed);
randn('state', seed);
ncases = 1000;
data = cell(N, ncases);
for m=1:ncases
  data(:,m) = sample_bnet(bnet);
end

priors = [0.1 1 10];
P = length(priors);
params = cell(1,P);
for p=1:P
  params{p} = cell(1,N);
  for i=1:N
    %params{p}{i} = {'prior', priors(p)};
    params{p}{i} = {'prior_type', 'dirichlet', 'dirichlet_weight', priors(p)};
  end
end

%sz = 1000:1000:10000;
sz = 10:10:100;
S = length(sz);
bic_score = zeros(S, 1);
bayes_score = zeros(S, P);
for i=1:S
  bic_score(i) = score_dags(data(:,1:sz(i)), ns, {dag}, 'scoring_fn', 'bic', 'params', []);
end
diff = zeros(S,P);
for p=1:P
  for i=1:S
    bayes_score(i,p) = score_dags(data(:,1:sz(i)), ns, {dag}, 'params', params{p});
  end
end

for p=1:P
  for i=1:S
    diff(i,p) = bayes_score(i,p)/ bic_score(i);
    %diff(i,p) = abs(bayes_score(i,p) - bic_score(i));
  end
end

if 0
plot(sz, diff(:,1), 'g--*', sz, diff(:,2), 'b-.+', sz, diff(:,3), 'k:s');
title('Relative BIC error vs. size of data set')
legend('BDeu 0.1', 'BDeu 1', 'Bdeu 10', 2)
end

if 0
plot(sz, bic_score, 'r-o',  sz, bayes_score(:,1), 'g--*', sz, bayes_score(:,2), 'b-.+', sz, bayes_score(:,3), 'k:s');
legend('bic', 'BDeu 0.01', 'BDeu 1', 'Bdeu 100')
ylabel('score')
title('score vs. size of data set')
end

%xlabel('num. data cases')

%previewfig(gcf, 'format', 'png', 'height', 2, 'color', 'rgb')
%exportfig(gcf, '/home/cs/murphyk/public_html/Bayes/Figures/bic.png', 'format', 'png', 'height', 2, 'color', 'rgb')
