% Make a QMR-like network 
% This is a bipartite graph, where the top layer contains hidden disease nodes,
% and the bottom later contains observed finding nodes.
% The diseases have Bernoulli CPDs, the findings noisy-or CPDs.
% See quickscore_inf_engine for references.

pMax = 0.01;
Nfindings = 10;
Ndiseases = 5;
%Nfindings = 20;
%Ndiseases = 10;

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

% Make the bnet in the straightforward way
tabular_leaves = 0;
obs_nodes = myunion(pos, neg) + Ndiseases;
big_bnet = mk_qmr_bnet(G, inhibit, leak, prior, tabular_leaves, obs_nodes);
big_evidence = cell(1, N);
big_evidence(findings(pos)) = num2cell(repmat(2, 1, length(pos)));
big_evidence(findings(neg)) = num2cell(repmat(1, 1, length(neg)));

%clf;draw_layout(big_bnet.dag);
%filename = '../public_html/Bayes/Figures/qmr.rnd.jpg';
%% 3x3 inches
%set(gcf,'units','inches');
%set(gcf,'PaperPosition',[0 0 3 3])  
%print(gcf,'-djpeg','-r100',filename);


% Marginalize out hidden leaves apriori
positive_leaves_only = 1;
[bnet, vals] = mk_minimal_qmr_bnet(G, inhibit, leak, prior, pos, neg, positive_leaves_only);
obs_nodes = bnet.observed;
evidence = cell(1, Ndiseases + length(obs_nodes));
evidence(obs_nodes) = num2cell(vals);


clear engine;
engine{1} = quickscore_inf_engine(inhibit, leak, prior);
engine{2} = jtree_inf_engine(big_bnet);
engine{3} = jtree_inf_engine(bnet);

%fname = '/home/cs/murphyk/matlab/Misc/loopybel.txt';
global BNT_HOME
fname = sprintf('%s/loopybel.txt', BNT_HOME);


max_iter = 6;
engine{4} = pearl_inf_engine(bnet, 'protocol', 'parallel', 'max_iter', max_iter);
%engine{5} = belprop_inf_engine(bnet, 'max_iter', max_iter, 'filename', fname);
engine{5} = belprop_inf_engine(bnet, 'max_iter', max_iter);

E = length(engine);
exact = 1:3;
loopy = [4 5];

ll = zeros(1,E);
tic; engine{1} = enter_evidence(engine{1}, pos, neg); toc
tic; [engine{2}, ll(2)] = enter_evidence(engine{2}, big_evidence); toc
tic; [engine{3}, ll(3)] = enter_evidence(engine{3}, evidence); toc
tic; [engine{4}, ll(4), niter(4)] = enter_evidence(engine{4}, evidence); toc
tic; [engine{5}, niter(5)] = enter_evidence(engine{5}, evidence); toc

ll

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

a = zeros(Ndiseases, 2);
for ei=1:length(loopy)
  for i=diseases(:)'
    a(i,ei) = approxeq(post(1, i), post(loopy(ei), i));
  end
end
disp('is the loopy posterior correct?');
disp(a)
