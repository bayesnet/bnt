% Bayesian model selection demo.

% We generate data from the model A->B
% and compute the posterior prob of all 3 dags on 2 nodes:
%  (1) A B,  (2) A <- B , (3) A -> B
% Models 2 and 3 are Markov equivalent, and therefore indistinguishable from 
% observational data alone.
% Using the "difficult" params, the true model only gets a higher posterior after 2000 trials!
% However, using the noisy NOT gate, the true model wins after 12 trials.

% ground truth
N = 2;
dag = zeros(N);
A = 1; B = 2; 
dag(A,B) = 1;

difficult = 0;
if difficult
  ntrials = 2000;
  ns = 3*ones(1,N);
  true_bnet = mk_bnet(dag, ns);
  rand('state', 0);
  temp = 5;
  for i=1:N
    %true_bnet.CPD{i} = tabular_CPD(true_bnet, i, temp);
    true_bnet.CPD{i} = tabular_CPD(true_bnet, i);
  end
else
  ntrials = 25;
  ns = 2*ones(1,N);
  true_bnet = mk_bnet(dag, ns);
  true_bnet.CPD{1} = tabular_CPD(true_bnet, 1, [0.5 0.5]);
  pfail = 0.1;
  psucc = 1-pfail;
  true_bnet.CPD{2} = tabular_CPD(true_bnet, 2, [pfail psucc; psucc pfail]); % NOT gate
end

G = mk_all_dags(N);
nhyp = length(G);
hyp_bnet = cell(1, nhyp);
for h=1:nhyp
  hyp_bnet{h} = mk_bnet(G{h}, ns);
  for i=1:N
    % We must set the CPTs to the mean of the prior for sequential log_marg_lik to be correct
    % The BDeu prior is score equivalent, so models 2,3 will be indistinguishable.
    % The uniform Dirichlet prior is not score equivalent...
    fam = family(G{h}, i);
    hyp_bnet{h}.CPD{i}= tabular_CPD(hyp_bnet{h}, i, 'prior_type', 'dirichlet', ...
				    'CPT', 'unif');
  end
end
prior = normalise(ones(1, nhyp));

% save results before doing sequential updating
init_hyp_bnet = hyp_bnet; 
init_prior = prior;


rand('state', 0);
hyp_w = zeros(ntrials+1, nhyp);
hyp_w(1,:) = prior(:)';

data = zeros(N, ntrials);

% First we compute the posteriors sequentially

LL = zeros(1, nhyp);
ll = zeros(1, nhyp);
for t=1:ntrials
  ev = cell2num(sample_bnet(true_bnet));
  data(:,t) = ev;
  for i=1:nhyp
    ll(i) = log_marg_lik_complete(hyp_bnet{i}, ev);
    hyp_bnet{i} = bayes_update_params(hyp_bnet{i}, ev);
  end
  prior = normalise(prior .* exp(ll));
  LL = LL + ll;
  hyp_w(t+1,:) = prior;
end

% Plot posterior model probabilities
% Red = model 1 (no arcs), blue/green = models 2/3 (1 arc)
% Blue = model 2 (2->1)
% Green = model 3 (1->2, "ground truth")

if 1
  figure;
m = size(hyp_w, 1);
h=plot(1:m, hyp_w(:,1), 'r-',  1:m, hyp_w(:,2), 'b-.', 1:m, hyp_w(:,3), 'g:');
axis([0 m   0 1])
title('model posterior vs. time')
%previewfig(gcf, 'format', 'png', 'height', 2, 'color', 'rgb')
%exportfig(gcf, '/home/cs/murphyk/public_html/Bayes/Figures/model_select.png',...
%'format', 'png', 'height', 2, 'color', 'rgb')
drawnow
end


% Now check that batch updating gives same result
hyp_bnet2 = init_hyp_bnet;
prior2 = init_prior;

cases = num2cell(data);
LL2 = zeros(1, nhyp);
for i=1:nhyp
  LL2(i) = log_marg_lik_complete(hyp_bnet2{i}, cases);
  hyp_bnet2{i} = bayes_update_params(hyp_bnet2{i}, cases);
end


assert(approxeq(LL, LL2))
LL

for i=1:nhyp
  for j=1:N
    s1 = struct(hyp_bnet{i}.CPD{j});
    s2 = struct(hyp_bnet2{i}.CPD{j});
    assert(approxeq(s1.CPT, s2.CPT))
  end
end

