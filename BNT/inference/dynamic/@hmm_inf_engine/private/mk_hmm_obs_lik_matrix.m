function obslik = mk_hmm_obs_lik_matrix(engine, evidence)

T  = size(evidence,2);
Q = length(engine.startprob);
obslik = ones(Q, T);
bnet = bnet_from_engine(engine);
% P(o1,o2| Q1,Q2) = P(o1|Q1,Q2) * P(o2|Q1,Q2)
onodes = bnet.observed;
for i=1:length(onodes)
  data = cell2num(evidence(onodes(i),:));
  if bnet.auto_regressive(onodes(i))
    params = engine.obsprob{i};
    mu = params.big_mu;
    Sigma = params.big_Sigma,
    W = params.big_W;
    mu0 = params.big_mu0;
    Sigma0 = params.big_Sigma0;
    %obslik_i = mk_arhmm_obs_lik(data, mu, Sigma, W, mu0, Sigma0
    obslik_i = clg_prob(data(:,1:T-1), data(:,2:T), mu, Sigma, W);
    obslik_i = [mixgauss_prob(data(:,1), mu0, Sigma0) obslik_i];
  elseif myismember(onodes(i), bnet.dnodes)
    %obslik_i = eval_pdf_cond_multinomial(data, engine.obsprob{i}.big_CPT);
    obslik_i = multinomial_prob(data, engine.obsprob{i}.big_CPT);
  else
    %obslik_i = eval_pdf_cond_gauss(data, engine.obsprob{i}.big_mu, engine.obsprob{i}.big_Sigma);
    obslik_i = mixgauss_prob(data, engine.obsprob{i}.big_mu, engine.obsprob{i}.big_Sigma);
  end
  obslik = obslik .* obslik_i;
end

