% influence diagram with no loops
%
% rv  dec
%  \  /
%  utility

N = 3;
dag = zeros(N);
X = 1; D = 2; U = 3;
dag([X D], U)=1;

ns = zeros(1,N);
ns(X) = 2; ns(D) = 2; ns(U) = 1;

limid = mk_limid(dag, ns, 'chance', X, 'decision', D, 'utility', U);

% use random params
limid.CPD{X} = tabular_CPD(limid, X);
limid.CPD{D} = tabular_decision_node(limid, D);
limid.CPD{U} = tabular_utility_node(limid, U);

%fname = '/home/cs/murphyk/matlab/Misc/loopybel.txt';
global BNT_HOME
fname = sprintf('%s/loopybel.txt', BNT_HOME);

engines = {};
engines{end+1} = global_joint_inf_engine(limid);
engines{end+1} = jtree_limid_inf_engine(limid);
%engines{end+1} = belprop_inf_engine(limid, 'max_iter', 2*N, 'filename', fname);
engines{end+1} = belprop_inf_engine(limid, 'max_iter', 2*N);

exact = [1 2];
approx = 3;

E = length(engines);
strategy = cell(1, E);
MEU = zeros(1, E);
for e=1:E
  [strategy{e}, MEU(e)] = solve_limid(engines{e});
  MEU
end
MEU

for e=exact(:)'
  assert(approxeq(strategy{exact(1)}{D}, strategy{e}{D}))
end

for e=approx(:)'
  approxeq(strategy{exact(1)}{D}, strategy{e}{D})
end
