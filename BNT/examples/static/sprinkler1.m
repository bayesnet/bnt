% Lawn sprinker example from Russell and Norvig p454
% For a picture, see http://www.cs.berkeley.edu/~murphyk/Bayes/usage.html#basics

N = 4; 
dag = zeros(N,N);
C = 1; S = 2; R = 3; W = 4;
dag(C,[R S]) = 1;
dag(R,W) = 1;
dag(S,W)=1;

false = 1; true = 2;
ns = 2*ones(1,N); % binary nodes

%bnet = mk_bnet(dag, ns);
bnet = mk_bnet(dag, ns, 'names', {'cloudy','S','R','W'}, 'discrete', 1:4);
names = bnet.names;
%C = names{'cloudy'};
bnet.CPD{C} = tabular_CPD(bnet, C, [0.5 0.5]);
bnet.CPD{R} = tabular_CPD(bnet, R, [0.8 0.2 0.2 0.8]);
bnet.CPD{S} = tabular_CPD(bnet, S, [0.5 0.9 0.5 0.1]);
bnet.CPD{W} = tabular_CPD(bnet, W, [1 0.1 0.1 0.01 0 0.9 0.9 0.99]);


CPD{C} = reshape([0.5 0.5], 2, 1);
CPD{R} = reshape([0.8 0.2 0.2 0.8], 2, 2);
CPD{S} = reshape([0.5 0.9 0.5 0.1], 2, 2);
CPD{W} = reshape([1 0.1 0.1 0.01 0 0.9 0.9 0.99], 2, 2, 2);
joint = zeros(2,2,2,2);
for c=1:2
  for r=1:2
    for s=1:2
      for w=1:2
	joint(c,s,r,w) = CPD{C}(c) * CPD{S}(c,s) * CPD{R}(c,r) * ...
	    CPD{W}(s,r,w);
      end
    end
  end
end

joint2 = repmat(reshape(CPD{C}, [2 1 1 1]), [1 2 2 2]) .* ...
	 repmat(reshape(CPD{S}, [2 2 1 1]), [1 1 2 2]) .* ...
	 repmat(reshape(CPD{R}, [2 1 2 1]), [1 2 1 2]) .* ...
	 repmat(reshape(CPD{W}, [1 2 2 2]), [2 1 1 1]);

assert(approxeq(joint, joint2));


engine = jtree_inf_engine(bnet);

evidence = cell(1,N);
evidence{W} = true;

[engine, ll] = enter_evidence(engine, evidence);

m = marginal_nodes(engine, S);
p1 = m.T(true) % P(S=true|W=true) = 0.4298
lik1 = exp(ll); % P(W=true) = 0.6471
assert(approxeq(p1, 0.4298));
assert(approxeq(lik1, 0.6471));

pSandW = sumv(joint(:,true,:,true), [C R]); % P(S,W) = sum_cr P(CSRW)
pW = sumv(joint(:,:,:,true), [C S R]);
pSgivenW = pSandW / pW; % P(S=t|W=t) = P(S=t,W=t)/P(W=t)
assert(approxeq(pW, lik1))
assert(approxeq(pSgivenW, p1))


m = marginal_nodes(engine, R);
p2 = m.T(true)  % P(R=true|W=true) =  0.7079     

pRandW = sumv(joint(:,:,true,true), [C S]); % P(R,W) = sum_cr P(CSRW)
pRgivenW = pRandW / pW; % P(R=t|W=t) = P(R=t,W=t)/P(W=t)
assert(approxeq(pRgivenW, p2))


% Add extra evidence that R=true
evidence{R} = true;

[engine, ll] = enter_evidence(engine, evidence);

m = marginal_nodes(engine, S);
p3 = m.T(true) % P(S=true|W=true,R=true) = 0.1945 
assert(approxeq(p3, 0.1945))


pSandRandW = sumv(joint(:,true,true,true), [C]); % P(S,R,W) = sum_c P(cSRW)
pRandW = sumv(joint(:,:,true,true), [C S]); % P(R,W) = sum_cs P(cSRW)
pSgivenWR = pSandRandW / pRandW; % P(S=t|W=t,R=t) = P(S=t,R=t,W=t)/P(W=t,R=t)
assert(approxeq(pSgivenWR, p3))

% So the sprinkler is less likely to be on if we know that
% it is raining, since the rain can "explain away" the fact
% that the grass is wet.

lik3 = exp(ll); % P(W=true, R=true) = 0.4581
% So the combined evidence is less likely (of course)




% Joint distributions

evidence = cell(1,N);
[engine, ll] = enter_evidence(engine, evidence);
m = marginal_nodes(engine, [S R W]);

evidence{R} = 2;
[engine, ll] = enter_evidence(engine, evidence);
m = marginal_nodes(engine, [S R W]);



