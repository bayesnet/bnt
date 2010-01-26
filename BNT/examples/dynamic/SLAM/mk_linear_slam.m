function [A,B,C,Q,R,Qbig,Rbig,init_x,init_V,robot_block,landmark_block,...
	  true_landmark_pos, true_robot_pos, true_data_assoc, ...
	  obs_rel_pos, ctrl_signal] = mk_linear_slam(varargin)

% We create data from a linear system for testing SLAM algorithms.
% i.e. , new robot pos = old robot pos + ctrl_signal, which is just a displacement vector.
% and  observation = landmark_pos - robot_pos, which is just a displacement vector.
%
% The behavior is determined by the following optional arguments:
%
% 'nlandmarks' - num. landmarks
% 'landmarks' - 'rnd' means random locations in the unit sqyare
%               'square' means at [1 1], [4 1], [4 4] and [1 4]
% 'T' - num steps to run
% 'ctrl' - 'stationary' means the robot remains at [0 0],
%          'leftright' means the robot receives a constant contol of [1 0],
%          'square' means we navigate the robot around the square
% 'data-assoc' - 'rnd' means we observe landmarks at random
%                'nn' means we observe the nearest neighbor landmark
%                'cycle' means we observe landmarks in order 1,2,.., 1, 2, ...

args = varargin;
% get mandatory params
for i=1:2:length(args)
  switch args{i},
   case 'nlandmarks', nlandmarks = args{i+1};
   case 'T', T = args{i+1};
  end
end

% set defaults
true_landmark_pos = rand(2,nlandmarks);
true_data_assoc = [];

% get args
for i=1:2:length(args)
  switch args{i},
   case 'landmarks',
    switch args{i+1},
     case 'rnd',   true_landmark_pos = rand(2,nlandmarks);
     case 'square',   true_landmark_pos = [1 1; 4 1; 4 4; 1 4]';
    end
   case 'ctrl',
    switch args{i+1},
     case 'stationary', ctrl_signal = repmat([0 0]', 1, T);
     case 'leftright', ctrl_signal = repmat([1 0]', 1, T);
     case 'square',   ctrl_signal = [repmat([1 0]', 1, T/4) repmat([0 1]', 1, T/4) ...
		    repmat([-1 0]', 1, T/4) repmat([0 -1]', 1, T/4)];
    end
   case 'data-assoc', 
    switch args{i+1},
     case 'rnd', true_data_assoc  = sample_discrete(normalise(ones(1,nlandmarks)),1,T);
     case 'cycle', true_data_assoc = wrap(1:T, nlandmarks);
    end
  end
end
if isempty(true_data_assoc)
  use_nn = 1;
else
  use_nn = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%
% generate data

init_robot_pos = [0 0]';
true_robot_pos = zeros(2, T);
true_rel_dist = zeros(2, T);
for t=1:T
  if t>1
    true_robot_pos(:,t) = true_robot_pos(:,t-1) + ctrl_signal(:,t);
  else
    true_robot_pos(:,t) = init_robot_pos + ctrl_signal(:,t);
  end
  nn = argmin(dist2(true_robot_pos(:,t)', true_landmark_pos'));
  if use_nn
    true_data_assoc(t) = nn;
  end
  true_rel_dist(:,t) = true_landmark_pos(:, nn) - true_robot_pos(:,t);
end


R = 1e-3*eye(2); % noise added to observation
Q = 1e-3*eye(2); % noise added to robot motion

% Create data set
obs_noise_seq = sample_gaussian([0 0]', R, T)';
obs_rel_pos = true_rel_dist + obs_noise_seq;
%obs_rel_pos = true_rel_dist;

%%%%%%%%%%%%%%%%%%
% Create params


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
%init_V(bi, bi) = 1e-5*eye(2); % very sure of robot posn
init_V(bi, bi) = Q; % simualate uncertainty due to 1 motion step
for i=1:nlandmarks
  bi = landmark_block(:,i);
  init_V(bi,bi)= 1e5*eye(2); % very uncertain of landmark psosns
  %init_x(bi) = true_landmark_pos(:,i);
  %init_V(bi,bi)= 1e-5*eye(2); % very sure of landmark psosns
end
