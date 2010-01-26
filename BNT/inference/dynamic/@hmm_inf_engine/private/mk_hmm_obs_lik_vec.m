function obslik = mk_hmm_obs_lik_vec(engine, evidence)

% P(o1,o2| h) = P(o1|h) * P(o2|h) where h = Q1,Q2,...

bnet = bnet_from_engine(engine);
ss = length(bnet.intra);
onodes = bnet.observed;
hnodes = mysetdiff(1:ss, onodes);
ns = bnet.node_sizes(:);
ns(onodes) = 1;

Q = length(engine.startprob);
obslik = ones(Q, 1);

for i=1:length(onodes)
  o = onodes(i);
  %data = cell2num(evidence(o,1));
  data = evidence{o,1};
  if myismember(o, bnet.dnodes)
    obslik_i = eval_pdf_cond_multinomial(data, engine.obsprob{i}.CPT);
  else
    if bnet.auto_regressive(o)
      error('can''t handle AR nodes')
    end
    %% calling mk_ghmm_obs_lik, which calls gaussian_prob, is slow, so we inline it
    %% and use the pre-computed  inverse matrix
    %obslik_i = mk_ghmm_obs_lik(data, engine.obsprob{i}.mu, engine.obsprob{i}.Sigma);
    x = data(:);
    m = engine.obsprob{i}.mu;
    Qi = size(m, 2);
    obslik_i = size(Qi, 1);
    invC = engine.obsprob{i}.inv_Sigma;
    denom = engine.obsprob{i}.denom;
    for j=1:Qi
      numer = exp(-0.5 * (x-m(:,j))' * invC(:,:,j) * (x-m(:,j)));
      obslik_i(j) = numer / denom(j);
    end
  end
  % convert P(o|ps) into P(o|h) by multiplying onto a (h,o) potential of all 1s
  ps = bnet.parents{o};
  dom = [ps o];
  obspot_i = dpot(dom, ns(dom), obslik_i);
  dom = [hnodes o];
  obspot = dpot(dom, ns(dom));
  obspot = multiply_by_pot(obspot, obspot_i);
  % compute p(oi|h) * p(oj|h)
  S = struct(obspot);
  obslik = obslik .* S.T(:);
end



