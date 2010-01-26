% This is like skf_data_assoc_gmux, except the objects don't move.
% We are uncertain of their initial positions, and get more and more observations
% over time. The goal is to test deterministic links (0 covariance).
% This is like robot1, except the robot doesn't move and is always at [0 0],
% so the relative location is simply L(s).

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

Xsz = 2; % state space = (x y)
Ysz = 2;
ns = zeros(1,N);
ns(Xs) = Xsz;
ns(Y) = Ysz;
ns(S) = nobj;

bnet = mk_dbn(intra, inter, ns, 'discrete', S, 'observed', [S Y]);

% For each object, we have
% X(t+1) = F X(t) + noise(Q)
% Y(t) = H X(t) + noise(R)
F = eye(2);
H = eye(2);
Q = 0*eye(Xsz); % no noise in dynamics
R = eye(Ysz);

init_state{1} = [10 10]';
init_state{2} = [10 -10]';
init_cov = eye(2);

% Uncertain of initial state (position)
for i=1:nobj
  bnet.CPD{Xs(i)} = gaussian_CPD(bnet, Xs(i), 'mean', init_state{i}, 'cov', init_cov);
end
bnet.CPD{S} = root_CPD(bnet, S); % always observed
bnet.CPD{Y} = gmux_CPD(bnet, Y, 'cov', repmat(R, [1 1 nobj]), 'weights', repmat(H, [1 1 nobj]));
% slice 2
eclass = bnet.equiv_class;
for i=1:nobj
  bnet.CPD{eclass(Xs(i), 2)} = gaussian_CPD(bnet, Xs(i)+N, 'mean', zeros(Xsz,1), 'cov', Q, 'weights', F);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create LDS params 

% X(t) = A X(t-1) + B U(t) + noise(Q)

% [L11]  = [1  ]  * [L1]       +   [Q ]
% [L2]     [  1]    [L2]           [ Q]

% Y(t)|S(t)=s  = C(s) X(t) + noise(R)
% Yt|St=1 = [1 0]  * [L1]  + R
%                    [L2]    

nlandmarks = nobj;

% Create indices into block structure
bs = 2*ones(1, nobj); % sizes of blocks in state space
for i=1:nlandmarks
  landmark_block(:,i) = block(i, bs)';
end
Xsz = 2*(nlandmarks); % 2 values for each landmark plus robot
Ysz = 2; % observe relative location

% create block-diagonal trans matrix for each switch
A = zeros(Xsz, Xsz);
for i=1:nlandmarks
  bi = landmark_block(:,i);
  A(bi, bi) = eye(2);
end
A = repmat(A, [1 1 nlandmarks]); % same for all switch values

% create block-diagonal system cov
Qbig = zeros(Xsz, Xsz);
Qbig = repmat(Qbig, [1 1 nlandmarks]);


% create observation matrix for each value of the switch node
% C(:,:,i) = (0 ... I ...) where the I is in the i'th posn.
C = zeros(Ysz, Xsz, nlandmarks);
for i=1:nlandmarks
  C(:, landmark_block(:,i), i) = eye(2);
end

% create observation cov for each value of the switch node
Rbig = repmat(R, [1 1 nlandmarks]);

% initial conditions
init_x = [init_state{1}; init_state{2}];
init_V = zeros(Xsz, Xsz);
for i=1:nlandmarks
  bi = landmark_block(:,i);
  init_V(bi,bi) = init_cov;
end



%%%%%%%%%%%%%%%%
% Observe objects at random
T = 10;
evidence = cell(N, T);
data_assoc = sample_discrete(normalise(ones(1,nobj)), 1, T);
evidence(S,:) = num2cell(data_assoc);
evidence = sample_dbn(bnet, 'evidence', evidence);


% Inference
ev = cell(N,T);
ev(bnet.observed,:) = evidence(bnet.observed, :);
y = cell2num(evidence(Y,:));

engine = pearl_unrolled_dbn_inf_engine(bnet);
engine = enter_evidence(engine, ev);

loopy_est_pos = zeros(2, nlandmarks);
loopy_est_pos_cov = zeros(2, 2, nlandmarks);
for i=1:nobj
  m = marginal_nodes(engine, Xs(i), T);
  loopy_est_pos(:,i) = m.mu;
  loopy_est_pos_cov(:,:,i) = m.Sigma;
end


[xsmooth, Vsmooth] = kalman_smoother(y, A, C, Qbig, Rbig, init_x, init_V, 'model', data_assoc);

kf_est_pos = zeros(2, nlandmarks);
kf_est_pos_cov = zeros(2, 2, nlandmarks);
for i=1:nlandmarks
  bi = landmark_block(:,i);
  kf_est_pos(:,i) = xsmooth(bi, T);
  kf_est_pos_cov(:,:,i) = Vsmooth(bi, bi, T);
end


kf_est_pos
loopy_est_pos

kf_est_pos_time = zeros(2, nlandmarks, T);
for t=1:T
  for i=1:nlandmarks
    bi = landmark_block(:,i);
    kf_est_pos_time(:,i,t) = xsmooth(bi, t);
  end
end
kf_est_pos_time % same for all t since smoothed
