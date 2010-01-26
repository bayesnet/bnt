function bnet = mk_square_hhmm(discrete_obs, true_params, topright)

% Make a 3 level  HHMM described by the following grammar
%
% Square -> CLK | CCK % clockwise or counterclockwise
% CLK -> LR UD RL DU start on top left (1 2 3 4)
% CCK -> RL UD LR DU  if start at top right (3 2 1 4)
% CCK -> UD LR DU RL if start at top left (2 1 4 3)
%
% LR = left-right, UD = up-down, RL = right-left, DU = down-up
% LR, UD, RL, DU are sub HMMs.
%
% For discrete observations, the subHMMs are 2-state left-right.
% LR emits L then l, etc.
%
% For cts observations, the subHMMs are 1 state.
% LR emits a vector in the -> direction, with a little noise.
% Since there is no constraint that we remain in the LR state as long as the RL state,
% the sides of the square might have different lengths,
% so the result is not really a square!
%
% If true_params = 0, we use random parameters at the top 2 levels
% (ready for learning). At the bottom level, we use noisy versions
% of the "true" observations.
%
% If topright=1, counter-clockwise starts at top right, not top left
% This example was inspired by Ivanov and Bobick.

if nargin < 3, topright = 1; end

if 1 % discrete_obs
  Qsizes = [2 4 2];
else
  Qsizes = [2 4 1];
end

D = 3;
Qnodes = 1:D;
startprob = cell(1,D);
transprob = cell(1,D);
termprob = cell(1,D);

% LEVEL 1

startprob{1} = 'unif';
transprob{1} = 'unif';

% LEVEL 2

if true_params
  startprob{2} = zeros(2, 4);
  startprob{2}(1, :) = [1 0 0 0];
  if topright
    startprob{2}(2, :) = [0 0 1 0];
  else
    startprob{2}(2, :) = [0 1 0 0];
  end
  
  transprob{2} = zeros(4, 2, 4);
  
  transprob{2}(:,1,:) = [0 1 0 0
		    0 0 1 0
		    0 0 0 1
		    0 0 0 1]; % 4->e
  if topright
    transprob{2}(:,2,:) = [0 0 0 1
		    1 0 0 0
		    0 1 0 0
		    0 0 0 1]; % 4->e
  else
    transprob{2}(:,2,:) = [0 0 0 1
		    1 0 0 0
		    0 0 1 0 % 3->e
		    0 0 1 0];
  end
  
  %termprob{2} = 'rightstop';
  termprob{2} = zeros(2,4);
  pfin = 0.8;
  termprob{2}(1,:) = [0 0 0 pfin]; % finish in state 4 (DU)
  if topright
    termprob{2}(2,:) = [0 0 0 pfin];
  else
    termprob{2}(2,:) = [0 0 pfin 0];  % finish in state 3 (RL)
  end
else
  % In the unsupervised case, it is essential that we break symmetry
  % in the initial param estimates.
  %startprob{2} = 'unif';
  %transprob{2} = 'unif';
  %termprob{2} = 'unif';
  startprob{2} = 'rnd';
  transprob{2} = 'rnd';
  termprob{2} = 'rnd';
end

% LEVEL 3

if 1 |  true_params
  startprob{3} = 'leftstart';
  transprob{3}  = 'leftright';
  termprob{3} = 'rightstop';
else
  % If we want to be able to run a base-level model backwards...
  startprob{3} = 'rnd';
  transprob{3}  = 'rnd';
  termprob{3} = 'rnd';
end
 

% OBS LEVEl

if discrete_obs
  % Initialise observations of lowest level primitives in a way which we can interpret
  chars = ['L', 'l', 'U', 'u', 'R', 'r', 'D', 'd'];
  L=find(chars=='L'); l=find(chars=='l');
  U=find(chars=='U'); u=find(chars=='u');
  R=find(chars=='R'); r=find(chars=='r');
  D=find(chars=='D'); d=find(chars=='d');
  Osize = length(chars);
  
  if true_params
    p = 1; % makes each state fully observed
  else
    p = 0.9;
  end
  
  obsprob = (1-p)*ones([4 2 Osize]);
  %       Q2 Q3 O
  obsprob(1, 1, L) =  p;
  obsprob(1, 2, l) =  p;
  obsprob(2, 1, U) =  p;
  obsprob(2, 2, u) =  p;
  obsprob(3, 1, R) =  p;
  obsprob(3, 2, r) =  p;
  obsprob(4, 1, D) =  p;
  obsprob(4, 2, d) =  p;
  obsprob = mk_stochastic(obsprob);
  Oargs = {'CPT', obsprob};
else
  % Initialise means of lowest level primitives in a way which we can interpret
  % These means are little vectors in the east, south, west, north directions.
  % (left-right=east, up-down=south, right-left=west, down-up=north)
  Osize = 2;
  mu = zeros(2, Qsizes(2), Qsizes(3));
  scale = 3;
  if true_params
    noise = 0;
  else
    noise = 0.5*scale;
  end
  for q3=1:Qsizes(3)
    mu(:, 1, q3) = scale*[1;0] + noise*rand(2,1);
  end
  for q3=1:Qsizes(3)
    mu(:, 2, q3) = scale*[0;-1] + noise*rand(2,1);
  end
  for q3=1:Qsizes(3)
    mu(:, 3, q3) = scale*[-1;0] + noise*rand(2,1);
  end
  for q3=1:Qsizes(3)
    mu(:, 4, q3) = scale*[0;1] + noise*rand(2,1);
  end
  Sigma = repmat(reshape(scale*eye(2), [2 2 1 1 ]), [1 1 Qsizes(2) Qsizes(3)]);
  Oargs = {'mean', mu, 'cov', Sigma, 'cov_type', 'diag'};
end

if discrete_obs
  selfprob = 0.5;
else
  selfprob = 0.95;
  % If less than this, it won't look like a square
  % because it doesn't spend enough time in each state
  % Unfortunately, the variance on durations (lengths of each side)
  % is very large
end
bnet = mk_hhmm('Qsizes', Qsizes, 'Osize', Osize', 'discrete_obs', discrete_obs, ...
	       'Oargs', Oargs, 'Ops', Qnodes(2:3), 'selfprob', selfprob, ...
	       'startprob', startprob, 'transprob', transprob, 'termprob', termprob);

