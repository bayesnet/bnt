% Burglar alarm example

N = 5;
dag = zeros(N,N);
E = 1; B = 2; R = 3; A = 4; C = 5;
dag(E,[R A]) = 1;
dag(B,A) = 1;
dag(A,C)=1;

% true = state 1, false = state 2
ns = 2*ones(1,N); % binary nodes
bnet = mk_bnet(dag, ns);

bnet.CPD{E} = tabular_CPD(bnet, E, [0.1 0.9]);
bnet.CPD{B} = tabular_CPD(bnet, B, [0.01 0.99]);
%bnet.CPD{R} = tabular_CPD(bnet, R, [0.65 0.00001 0.35 0.99999]);
bnet.CPD{R} = tabular_CPD(bnet, R, [0.65 0.01 0.35 0.99]);
bnet.CPD{A} = tabular_CPD(bnet, A, [0.95 0.8 0.3 0.001 0.05 0.2 0.7 0.999]);
bnet.CPD{C} = tabular_CPD(bnet, C, [0.7 0.05 0.3 0.95]);


engine = jtree_inf_engine(bnet);
ev  = cell(1,N);
ev{C} = 1;
engine = enter_evidence(engine, ev);
mE = marginal_nodes(engine, E);
mB = marginal_nodes(engine, B);
fprintf('P(E|c)=%5.3f, P(B|c)=%5.3f\n', mE.T(1), mB.T(1))

ev{C} = 1;
ev{R} = 1;
engine = enter_evidence(engine, ev);
mE = marginal_nodes(engine, E);
mB = marginal_nodes(engine, B);
fprintf('P(E|c,r)=%5.3f, P(B|c,r)=%5.3f\n', mE.T(1), mB.T(1))


if 0
nsamples = 100;
samples = zeros(nsamples, 5);
for i=1:nsamples
  samples(i,:) = cell2num(sample_bnet(bnet))';
end
end
