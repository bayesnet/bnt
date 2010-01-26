% Conditional Gaussian network with vector-valued nodes and random params

ns = 2*ones(1,9);
bnet  = mk_incinerator_bnet(ns);

engines = {};
%engines{end+1} = stab_cond_gauss_inf_engine(bnet);
engines{end+1} = jtree_inf_engine(bnet);
engines{end+1} = cond_gauss_inf_engine(bnet);

[err, time] = cmp_inference_static(bnet, engines, 'singletons_only', 1);
