% We navigate a robot around a square using a fixed control policy and no noise.
% We assume the robot observes the relative distance to the nearest landmark.
% Everything is linear-Gaussian.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create toy data set

seed = 0;
rand('state', seed);
randn('state', seed);

if 1
  T = 20;
  ctrl_signal = [repmat([1 0]', 1, T/4) repmat([0 1]', 1, T/4) ...
		 repmat([-1 0]', 1, T/4) repmat([0 -1]', 1, T/4)];
else
  T = 5;
  ctrl_signal = repmat([1 0]', 1, T);
end

nlandmarks = 4;
true_landmark_pos = [1 1;
		     4 1;
		     4 4;
		     1 4]';
init_robot_pos = [0 0]';

true_robot_pos = zeros(2, T);
true_data_assoc = zeros(1, T);
true_rel_dist = zeros(2, T);
for t=1:T
  if t>1
    true_robot_pos(:,t) = true_robot_pos(:,t-1) + ctrl_signal(:,t);
  else
    true_robot_pos(:,t) = init_robot_pos + ctrl_signal(:,t);
  end
  nn = argmin(dist2(true_robot_pos(:,t)', true_landmark_pos'));
  %nn = t; % observe 1, 2, 3
  true_data_assoc(t) = nn;
  true_rel_dist(:,t) = true_landmark_pos(:, nn) - true_robot_pos(:,t);
end

figure(1);
%clf; 
hold on
%plot(true_landmark_pos(1,:), true_landmark_pos(2,:), '*');
for i=1:nlandmarks
  text(true_landmark_pos(1,i), true_landmark_pos(2,i), sprintf('L%d',i));
end
for t=1:T
  text(true_robot_pos(1,t), true_robot_pos(2,t), sprintf('%d',t));
end
hold off
axis([-1 6 -1 6])

R = 1e-3*eye(2); % noise added to observation
Q = 1e-3*eye(2); % noise added to robot motion

% Create data set
obs_noise_seq = sample_gaussian([0 0]', R, T)';
obs_rel_pos = true_rel_dist + obs_noise_seq;
%obs_rel_pos = true_rel_dist;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create params for inference

% X(t) = A X(t-1) + B U(t) + noise(Q)

% [L1]  = [1     ]  * [L1]       + [0]  * Ut  + [0   ]
% [L2]    [  1   ]    [L2]         [0]          [ 0  ]
% [R ]t   [     1]    [R ]t-1      [1]          [   Q]

% Y(t)|S(t)=s  = C(s) X(t) + noise(R)
% Yt|St=1 = [1 0 -1]  * [L1]  + R
%                       [L2]    
%                       [R ]    

% Create indices into block structure
bs = 2*ones(1, nlandmarks+1); % sizes of blocks in state space
robot_block =  block(nlandmarks+1, bs);
for i=1:nlandmarks
  landmark_block(:,i) = block(i, bs)';
end
Xsz = 2*(nlandmarks+1); % 2 values for each landmark plus robot
Ysz = 2; % observe relative location
Usz = 2; % input is (dx, dy)


% create block-diagonal trans matrix for each switch
A = zeros(Xsz, Xsz);
for i=1:nlandmarks
  bi = landmark_block(:,i);
  A(bi, bi) = eye(2);
end
bi = robot_block;
A(bi, bi) = eye(2);
A = repmat(A, [1 1 nlandmarks]); % same for all switch values

% create block-diagonal system cov


Qbig = zeros(Xsz, Xsz);
bi = robot_block;
Qbig(bi,bi) = Q; % only add noise to robot motion
Qbig = repmat(Qbig, [1 1 nlandmarks]);

% create input matrix
B = zeros(Xsz, Usz);
B(robot_block,:) = eye(2); % only add input to robot position
B = repmat(B, [1 1 nlandmarks]);

% create observation matrix for each value of the switch node
% C(:,:,i) = (0 ... I ... -I) where the I is in the i'th posn.
% This computes L(i) - R
C = zeros(Ysz, Xsz, nlandmarks);
for i=1:nlandmarks
  C(:, landmark_block(:,i), i) = eye(2); 
  C(:, robot_block, i) = -eye(2);
end

% create observation cov for each value of the switch node
Rbig = repmat(R, [1 1 nlandmarks]);

% initial conditions
init_x = zeros(Xsz, 1);
init_v = zeros(Xsz, Xsz);
bi = robot_block;
init_x(bi) = init_robot_pos;
init_V(bi, bi) = 1e-5*eye(2); % very sure of robot posn
for i=1:nlandmarks
  bi = landmark_block(:,i);
  init_V(bi,bi)= 1e5*eye(2); % very uncertain of landmark psosns
  %init_x(bi) = true_landmark_pos(:,i);
  %init_V(bi,bi)= 1e-5*eye(2); % very sure of landmark psosns
end

%%%%%%%%%%%%%%%%%%%%%
% Inference
if 1
[xsmooth, Vsmooth] = kalman_smoother(obs_rel_pos, A, C, Qbig, Rbig, init_x, init_V, ...
				     'model', true_data_assoc, 'u', ctrl_signal, 'B', B);

est_robot_pos = xsmooth(robot_block, :);
est_robot_pos_cov = Vsmooth(robot_block, robot_block, :);

for i=1:nlandmarks
  bi = landmark_block(:,i);
  est_landmark_pos(:,i) = xsmooth(bi, T);
  est_landmark_pos_cov(:,:,i) = Vsmooth(bi, bi, T);
end
end


if 0
figure(1); hold on
for i=1:nlandmarks
  h=plotgauss2d(est_landmark_pos(:,i), est_landmark_pos_cov(:,:,i));
  set(h, 'color', 'r')
end
hold off

hold on
for t=1:T
  h=plotgauss2d(est_robot_pos(:,t), est_robot_pos_cov(:,:,t));
  set(h,'color','r')
  h=text(est_robot_pos(1,t), est_robot_pos(2,2), sprintf('R%d', t));
  set(h,'color','r')
end
hold off
end


if 0
figure(3)
if 0
  for t=1:T
    imagesc(inv(Vsmooth(:,:,t)))
    colorbar
    fprintf('t=%d; press key to continue\n', t);
    pause
  end
else
  for t=1:T
    subplot(5,4,t)
    imagesc(inv(Vsmooth(:,:,t)))
  end
end
end





%%%%%%%%%%%%%%%%%
% DBN inference

if 1
  [bnet, Unode, Snode, Lnodes, Rnode, Ynode, Lsnode] = ...
      mk_gmux_robot_dbn(nlandmarks, Q, R, init_x, init_V, robot_block, landmark_block);
  engine = pearl_unrolled_dbn_inf_engine(bnet, 'max_iter', 50, 'filename', ...
					 '/home/eecs/murphyk/matlab/loopyslam.txt');
else
  [bnet, Unode, Snode, Lnodes, Rnode, Ynode] = ...
      mk_gmux2_robot_dbn(nlandmarks, Q, R, init_x, init_V, robot_block, landmark_block);
  engine = jtree_dbn_inf_engine(bnet);
end

nnodes = bnet.nnodes_per_slice;
evidence = cell(nnodes, T);
evidence(Ynode, :) = num2cell(obs_rel_pos, 1);
evidence(Unode, :) = num2cell(ctrl_signal, 1);
evidence(Snode, :) = num2cell(true_data_assoc);


[engine, ll, niter] = enter_evidence(engine, evidence);
niter

loopy_est_robot_pos = zeros(2, T);
for t=1:T
  m = marginal_nodes(engine, Rnode, t);
  loopy_est_robot_pos(:,t) = m.mu;
end

for i=1:nlandmarks
  m = marginal_nodes(engine, Lnodes(i), T);
  loopy_est_landmark_pos(:,i) = m.mu;
  loopy_est_landmark_pos_cov(:,:,i) = m.Sigma;
end


