
seed = 0;
rand('state', seed);
randn('state', seed);

discrete_obs = 1;
topright = 0;

Qsizes = [2 4 2];
D = 3;
Qnodes = 1:D;
startprob = cell(1,D);
transprob = cell(1,D);
termprob = cell(1,D);

% LEVEL 1

startprob{1} = 'ergodic';
transprob{1} = 'ergodic';

% LEVEL 2

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
termprob{2} = zeros(2,4,2);
pfin = 0.8;
termprob{2}(1,:,2) = [0 0 0 pfin]; % finish in state 4 (DU)
termprob{2}(1,:,1) = 1 - [0 0 0 pfin];
if topright
  termprob{2}(2,:,2) = [0 0 0 pfin];
  termprob{2}(2,:,1) = 1 - [0 0 0 pfin];
else
  termprob{2}(2,:,2) = [0 0 pfin 0];  % finish in state 3 (RL)
  termprob{2}(2,:,1) = 1 - [0 0 pfin 0];
end

% LEVEL 3

startprob{3} = 'leftstart';
transprob{3}  = 'leftright';
termprob{3} = 'rightstop';


% OBS LEVEl

if discrete_obs
  chars = ['L', 'l', 'U', 'u', 'R', 'r', 'D', 'd'];
  L=find(chars=='L'); l=find(chars=='l');
  U=find(chars=='U'); u=find(chars=='u');
  R=find(chars=='R'); r=find(chars=='r');
  D=find(chars=='D'); d=find(chars=='d');
  Osize = length(chars);
  
  obsprob = zeros([4 2 Osize]);
  %       Q2 Q3 O
  obsprob(1, 1, L) =  1.0;
  obsprob(1, 2, l) =  1.0;
  obsprob(2, 1, U) =  1.0;
  obsprob(2, 2, u) =  1.0;
  obsprob(3, 1, R) =  1.0;
  obsprob(3, 2, r) =  1.0;
  obsprob(4, 1, D) =  1.0;
  obsprob(4, 2, d) =  1.0;
  
  Oargs = {'CPT', obsprob};
else
  Osize = 2;
  mu = zeros(2, 4, 2);
  noise = 0;
  scale = 10;
  for q3=1:2
    mu(:, 1, q3) = scale*[1;0] + noise*rand(2,1);
  end
  for q3=1:2
    mu(:, 2, q3) = scale*[0;-1] + noise*rand(2,1);
  end
  for q3=1:2
    mu(:, 3, q3) = scale*[-1;0] + noise*rand(2,1);
  end
  for q3=1:2
    mu(:, 4, q3) = scale*[0;1] + noise*rand(2,1);
  end
  Sigma = repmat(reshape(0.01*eye(2), [2 2 1 1 ]), [1 1 4 2]);
  Oargs = {'mean', mu, 'cov', Sigma};
end

bnet = mk_hhmm('Qsizes', Qsizes, 'Osize', Osize', 'discrete_obs', discrete_obs, ...
	       'Oargs', Oargs, 'Ops', Qnodes(2:3), ...
	       'startprob', startprob, 'transprob', transprob, 'termprob', termprob);

if discrete_obs
  Tmax = 30;
else
  Tmax = 200;
end
usecell = ~discrete_obs;
Q1 = 1; Q2 = 2; Q3 = 3; F3 = 4; F2 = 5; Onode = 6;
Qnodes = [Q1 Q2 Q3]; Fnodes = [F2 F3];

for seqi=1:3
  evidence = sample_dbn(bnet, Tmax, usecell, 'stop_sampling_F2');      
  T = size(evidence, 2)
  if discrete_obs
    pretty_print_hhmm_parse(evidence, Qnodes, Fnodes, Onode, chars);
  else
    pos = zeros(2,T+1);
    delta = cell2num(evidence(Onode,:));
    clf
    hold on
    cols = {'r', 'g', 'k', 'b'};
    boundary = cell2num(evidence(F3,:))-1;
    coli = 1;
    for t=2:T+1
      pos(:,t) = pos(:,t-1) + delta(:,t-1);
      plot(pos(1,t), pos(2,t), sprintf('%c.', cols{coli}));
      if boundary(t-1)
	coli = coli + 1;
	coli = mod(coli-1, length(cols)) + 1;
      end
    end
    %plot(pos(1,:), pos(2,:), '.')
    %pretty_print_hhmm_parse(evidence, Qnodes, Fnodes, Onode, []);
    pause
  end
end

eclass = bnet.equiv_class;
S=struct(bnet.CPD{eclass(Q2,2)});







