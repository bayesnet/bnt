% Compare different loopy belief propagation algorithms on a graph with many loops

bnet = mk_asia_bnet('orig');

engines = {};
engines{end+1} = jtree_inf_engine(bnet);
engines{end+1} = pearl_inf_engine(bnet, 'protocol', 'parallel');
engines{end+1} = belprop_fg_inf_engine(bnet_to_fgraph(bnet));
engines{end+1} = belprop_inf_engine(bnet, 'protocol', 'parallel');

[time, engines] = cmp_inference_static(bnet, engines, 'maximize', 0, 'exact', 1, 'observed', [1 3 5], ...
				   'check_ll', 0, 'singletons_only', 1, 'check_converged', 2:4);

