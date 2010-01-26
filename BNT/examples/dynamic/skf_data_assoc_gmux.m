% We consider a switching Kalman filter of the kind studied
% by Zoubin Ghahramani, i.e., where the switch node determines
% which of the hidden chains we get to observe (data association).
% e.g., for n=2 chains
% 
% X1 -> X1
% | X2 -> X2
% \ |
%  v
%  Y
%  ^
%  |
%  S
%
% Y is a gmux (multiplexer) node, where S switches in one of the parents.
% We differ from Zoubin by not connecting the S nodes over time (which
% doesn't make sense for data association).
% Indeed, we assume the S nodes are always observed.
% 
%
% We will track 2 objects (points) moving in the plane, as in BNT/Kalman/tracking_demo.
% We will alternate between observing them.

nobj = 2;
N = nobj+2;
Xs = 1:nobj;
S = nobj+1;
Y = nobj+2;

intra = zeros(N,N);
inter = zeros(N,N);
intra([Xs S], Y) =1;
for i=1:nobj
  inter(Xs(i), Xs(i))=1;
end

Xsz = 4; % state space = (x y xdot ydot)
Ysz = 2;
ns = zeros(1,N);
ns(Xs) = Xsz;
ns(Y) = Ysz;
ns(S) = n;

bnet = mk_dbn(intra, inter, ns, 'discrete', S, 'observed', [S Y]);

% For each object, we have
% X(t+1) = F X(t) + noise(Q)
% Y(t) = H X(t) + noise(R)
F = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1];
H = [1 0 0 0; 0 1 0 0];
Q = 1e-3*eye(Xsz);
%R = 1e-3*eye(Ysz);
R = eye(Ysz);

% We initialise object 1 moving to the right, and object 2 moving to the left
% (Here, we assume nobj=2)
init_state{1} = [10 10 1 0]';
init_state{2} = [10 -10 -1 0]';

for i=1:nobj
  bnet.CPD{Xs(i)} = gaussian_CPD(bnet, Xs(i), 'mean', init_state{i}, 'cov', 1e-4*eye(Xsz));
end
bnet.CPD{S} = root_CPD(bnet, S); % always observed
bnet.CPD{Y} = gmux_CPD(bnet, Y, 'cov', repmat(R, [1 1 nobj]), 'weights', repmat(H, [1 1 nobj]));
% slice 2
eclass = bnet.equiv_class;
for i=1:nobj
  bnet.CPD{eclass(Xs(i), 2)} = gaussian_CPD(bnet, Xs(i)+N, 'mean', zeros(Xsz,1), 'cov', Q, 'weights', F);
end

% Observe objects at random
T = 10;
evidence = cell(N, T);
data_assoc = sample_discrete(normalise(ones(1,nobj)), 1, T);
evidence(S,:) = num2cell(data_assoc);
evidence = sample_dbn(bnet, 'evidence', evidence);

% plot the data
true_state = cell(1,nobj);
for i=1:nobj
  true_state{i} = cell2num(evidence(Xs(i), :)); % true_state{i}(:,t) = [x y xdot ydot]'
end
obs_pos = cell2num(evidence(Y,:));
figure(1)
clf
hold on
styles = {'rx', 'go', 'b+', 'k*'};
for i=1:nobj
  plot(true_state{i}(1,:), true_state{i}(2,:), styles{i});
end
for t=1:T
  text(obs_pos(1,t), obs_pos(2,t), sprintf('%d', t));
end
hold off
relax_axes(0.1)


% Inference
ev = cell(N,T);
ev(bnet.observed,:) = evidence(bnet.observed, :);

engines = {};
engines{end+1} = jtree_dbn_inf_engine(bnet);
%engines{end+1} = scg_unrolled_dbn_inf_engine(bnet, T);
engines{end+1} = pearl_unrolled_dbn_inf_engine(bnet);
E = length(engines);

inferred_state = cell(nobj,E); % inferred_state{i,e}(:,t)
for e=1:E
  engines{e} = enter_evidence(engines{e}, ev);
  for i=1:nobj
    inferred_state{i,e} = zeros(4, T);
    for t=1:T
      m = marginal_nodes(engines{e}, Xs(i), t);
      inferred_state{i,e}(:,t) = m.mu;
    end
  end
end
inferred_state{1,1}
inferred_state{1,2}

% Plot results
figure(2)
clf
hold on
styles = {'rx', 'go', 'b+', 'k*'};
nstyles = length(styles);
c = 1;
for e=1:E
  for i=1:nobj
    plot(inferred_state{i,e}(1,:), inferred_state{i,e}(2,:), styles{mod(c-1,nstyles)+1});
    c = c + 1;
  end
end
for t=1:T
  text(obs_pos(1,t), obs_pos(2,t), sprintf('%d', t));
end
hold off
relax_axes(0.1)
