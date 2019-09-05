if 0
% Generate some sample paths

bnet = mk_map_hhmm('p', 1);
% assign numbers to the nodes in topological order
U = 1; A = 2; C = 3; F = 4; O = 5;


seed = 0;
rand('state', seed);
randn('state', seed);

% control policy = sweep right then left
T = 10;
ss = 5;
ev = cell(ss, T);
ev(U,:) = num2cell([R*ones(1,5) L*ones(1,5)]);

% fix initial conditions to be in left most state
ev{A,1} = 1; 
ev{C,1} = 1; 
evidence = sample_dbn(bnet, 'length', T, 'evidence', ev)


% Now do same but with noisy actuators

bnet = mk_map_hhmm('p', 0.8);
evidence = sample_dbn(bnet, 'length', T, 'evidence', ev)

end

% Now do same but with 4 observations per slice

bnet = mk_map_hhmm('p', 0.8, 'obs_model', 'four');
ss = bnet.nnodes_per_slice;

ev = cell(ss, T);
ev(U,:) = num2cell([R*ones(1,5) L*ones(1,5)]);
ev{A,1} = 1; 
ev{C,1} = 1; 
evidence = sample_dbn(bnet, 'length', T, 'evidence', ev)
