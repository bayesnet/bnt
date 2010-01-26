function engine = var_elim_inf_engine(bnet, varargin)
% VAR_ELIM_INF_ENGINE Variable elimination inference engine
% engine = var_elim_inf_engine(bnet)
%
% For details on variable elimination, see
% - R. Dechter, "Bucket Elimination: A Unifying Framework for Probabilistic Inference", UA1 96, pp. 211-219. 
% - Z. Li and B. D'Ambrosio, "Efficient inference in Bayes networks as a combinatorial
%     optimization problem", Intl. J. Approximate Reasoning, 11(1):55-81, 1994
% - R. McEliece and S. M. Aji, "The Generalized Distributive Law", IEEE Trans. Inform. Theory, 46(2), 2000


% This is where we will store the results between enter_evidence and marginal_nodes
engine.evidence = [];

engine = class(engine, 'var_elim_inf_engine', inf_engine(bnet));
