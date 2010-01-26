function engine = gaussian_inf_engine(bnet)
% GAUSSIAN_INF_ENGINE Computes the joint multivariate Gaussian corresponding to the bnet
% engine = gaussian_inf_engine(bnet)
%
% For details on how to compute the joint Gaussian from the bnet, see
% - "Gaussian Influence Diagrams", R. Shachter and C. R. Kenley, Management Science, 35(5):527--550, 1989.
% Once we have the Gaussian, we can apply the standard formulas for conditioning and marginalization.

assert(isequal(bnet.cnodes, 1:length(bnet.dag)));

[W, D, mu] = extract_params_from_gbn(bnet);
U = inv(eye(size(W)) - W')';
Sigma = U' * D * U;

engine.mu = mu;
engine.Sigma = Sigma;
%engine.logp = log(normal_coef(Sigma));

% This is where we will store the results between enter_evidence and marginal_nodes  
engine.Hmu = [];
engine.HSigma = [];
engine.hnodes = [];

engine = class(engine, 'gaussian_inf_engine', inf_engine(bnet));

