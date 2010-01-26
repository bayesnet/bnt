% Online Bayesian model selection demo.

% We generate data from the model A->B
% and compute the posterior prob of all 3 dags on 2 nodes:
%  (1) A B,  (2) A <- B , (3) A -> B
% Models 2 and 3 are Markov equivalent, and therefore indistinguishable from 
% observational data alone.

% We control the dependence of B on A by setting
% P(B|A) = 0.5 - epislon and vary epsilon
% as in Koller & Friedman book p512

% ground truth
N = 2;
dag = zeros(N);
A = 1; B = 2; 
dag(A,B) = 1;

ntrials = 100;
ns = 2*ones(1,N);
true_bnet = mk_bnet(dag, ns);
true_bnet.CPD{1} = tabular_CPD(true_bnet, 1, [0.5 0.5]);

% hypothesis space
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

clf
seeds = 1:3;
expt = 1;
for seedi=1:length(seeds)
  seed = seeds(seedi);
  rand('state', seed);
  randn('state', seed);
    
  es = [0.05 0.1 0.15 0.2];
  for ei=1:length(es)
    e = es(ei);
    true_bnet.CPD{2} = tabular_CPD(true_bnet, 2, [0.5+e 0.5-e; 0.5-e 0.5+e]);

    prior = normalise(ones(1, nhyp));
    hyp_w = zeros(ntrials+1, nhyp);
    hyp_w(1,:) = prior(:)';
    LL = zeros(1, nhyp);
    ll = zeros(1, nhyp);
    for t=1:ntrials
      ev = cell2num(sample_bnet(true_bnet));
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
    
    subplot2(length(seeds), length(es), seedi, ei);
    m = size(hyp_w,1);
    h=plot(1:m, hyp_w(:,1), 'r-',  1:m, hyp_w(:,2), 'b-.', 1:m, hyp_w(:,3), 'g:');
    axis([0 m   0 1])
    %title('model posterior vs. time')
    title(sprintf('e=%3.2f, seed=%d', e, seed));
    drawnow
    expt = expt + 1;
  end
end
