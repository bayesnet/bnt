% Test jtree_compiled on a toy QMR network.

rand('state', 0);
randn('state', 0);
pMax = 0.01;
Nfindings = 10;
Ndiseases = 5;

N=Nfindings+Ndiseases;
findings = Ndiseases+1:N;
diseases = 1:Ndiseases;

G = zeros(Ndiseases, Nfindings);
for i=1:Nfindings
  v= rand(1,Ndiseases);
  rents = find(v<0.8);
  if (length(rents)==0)
    rents=ceil(rand(1)*Ndiseases);
  end
  G(rents,i)=1;
end       

prior = pMax*rand(1,Ndiseases);
leak = 0.5*rand(1,Nfindings); % in real QMR, leak approx exp(-0.02) = 0.98     
%leak = ones(1,Nfindings); % turns off leaks, which makes inference much harder
inhibit = rand(Ndiseases, Nfindings);
inhibit(not(G)) = 1;

% first half of findings are +ve, second half -ve
% The very first and last findings are hidden
pos = 2:floor(Nfindings/2);
neg = (pos(end)+1):(Nfindings-1);

big = 1;

if big
  % Make the bnet in the straightforward way
  tabular_leaves = 1;
  obs_nodes = myunion(pos, neg) + Ndiseases;
  bnet = mk_qmr_bnet(G, inhibit, leak, prior, tabular_leaves, obs_nodes);
  evidence = cell(1, N);
  evidence(findings(pos)) = num2cell(repmat(2, 1, length(pos)));
  evidence(findings(neg)) = num2cell(repmat(1, 1, length(neg)));
else
  % Marginalize out hidden leaves apriori
  positive_leaves_only = 1;
  [bnet, vals] = mk_minimal_qmr_bnet(G, inhibit, leak, prior, pos, neg, positive_leaves_only);
  obs_nodes = bnet.observed;
  evidence = cell(1, Ndiseases + length(obs_nodes));
  evidence(obs_nodes) = num2cell(vals);
end

engine = {};
engine{end+1} = jtree_inf_engine(bnet);

E = length(engine);
exact = 1:E;
ll = zeros(1,E);
for e=1:E
  tic; [engine{e}, ll(e)] = enter_evidence(engine{e}, evidence); toc
end

assert(all(approxeq(ll(exact), ll(exact(1)))))

post = zeros(E, Ndiseases);
for e=1:E
  for i=diseases(:)'
    m = marginal_nodes(engine{e}, i);
    post(e, i) = m.T(2);
  end
end
for e=exact(:)'
  for i=diseases(:)'
    assert(approxeq(post(1, i), post(e, i)));
  end
end

