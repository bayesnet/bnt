% Same as cg1, except we assume all discretes are observed,
% and use loopy for approximate inference.

ns = 2*ones(1,9);
F = 1; W = 2; E = 3; B = 4; C = 5; D = 6; Min = 7; Mout = 8; L = 9;
n = 9;
dnodes = [B F W];
cnodes = mysetdiff(1:n, dnodes);

%bnet  = mk_incinerator_bnet(ns);
bnet  = mk_incinerator_bnet;

bnet.observed = [dnodes E];

engines = {};
engines{end+1} = jtree_inf_engine(bnet);
engines{end+1} = pearl_inf_engine(bnet, 'protocol', 'parallel');
nengines = length(engines);


[time, engines] = cmp_inference_static(bnet, engines, 'maximize', 0, 'check_ll', 0, ...
				      'singletons_only', 0, 'exact', 1, 'check_converged', 2);
