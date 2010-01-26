% This is like robot1, except we only use a Kalman filter.
% The goal is to study how the precision matrix changes.

seed = 0;
rand('state', seed);
randn('state', seed);

if 0
  T = 20;
  ctrl_signal = [repmat([1 0]', 1, T/4) repmat([0 1]', 1, T/4) ...
		 repmat([-1 0]', 1, T/4) repmat([0 -1]', 1, T/4)];
else
  T = 12;
  ctrl_signal = repmat([1 0]', 1, T);
end

nlandmarks = 6;
if 0
  true_landmark_pos = [1 1;
		    4 1;
		    4 4;
		    1 4]';
else
  true_landmark_pos = 10*rand(2,nlandmarks);
end
figure(1); clf
hold on
for i=1:nlandmarks
  %text(true_landmark_pos(1,i), true_landmark_pos(2,i), sprintf('L%d',i));
  plot(true_landmark_pos(1,i), true_landmark_pos(2,i), '*')
end
hold off

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
  %nn = argmin(dist2(true_robot_pos(:,t)', true_landmark_pos'));
  nn = wrap(t, nlandmarks); % observe 1, 2, 3, 4, 1, 2, ...
  true_data_assoc(t) = nn;
  true_rel_dist(:,t) = true_landmark_pos(:, nn) - true_robot_pos(:,t);
end

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

[xsmooth, Vsmooth] = kalman_filter(obs_rel_pos, A, C, Qbig, Rbig, init_x, init_V, ...
				     'model', true_data_assoc, 'u', ctrl_signal, 'B', B);

est_robot_pos = xsmooth(robot_block, :);
est_robot_pos_cov = Vsmooth(robot_block, robot_block, :);

for i=1:nlandmarks
  bi = landmark_block(:,i);
  est_landmark_pos(:,i) = xsmooth(bi, T);
  est_landmark_pos_cov(:,:,i) = Vsmooth(bi, bi, T);
end



P = zeros(size(Vsmooth));
for t=1:T
  P(:,:,t) = inv(Vsmooth(:,:,t));
end

figure(1)
for t=1:T
  subplot(T/2,2,t)
  imagesc(P(1:2:end,1:2:end, t))
  colorbar
end

figure(2)
for t=1:T
  subplot(T/2,2,t)
  imagesc(Vsmooth(1:2:end,1:2:end, t))
  colorbar
end



% marginalize out robot position and then check structure
bi = landmark_block(:);
V = Vsmooth(bi,bi,T); 
P = inv(V);
P(1:2:end,1:2:end)
