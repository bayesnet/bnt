% Make the following HHMM
%
%    S1 <----------------------> S2
%    |                           |
%    |                           |
%   M1 -> M2 -> M3 -> end        B1 -> end
%
% where Mi represents the i'th letter in the motif
% and B is the background state.
% Si chooses between running the motif or the background.
% The Si and B states have self loops (not shown).

if 0
seed = 0;
rand('state', seed);
randn('state', seed);
end

chars = ['a', 'c', 'g', 't'];
Osize = length(chars);

motif_pattern = 'acca';
motif_length = length(motif_pattern);
Qsize = [2 motif_length];
Qnodes = 1:2;
D = 2;
transprob = cell(1,D);
termprob = cell(1,D);
startprob = cell(1,D);

% startprob{d}(k,j), startprob{1}(1,j)
% transprob{d}(i,k,j), transprob{1}(i,j)
% termprob{d}(k,j)


% LEVEL 1

startprob{1} = zeros(1, 2);
startprob{1} = [1 0]; % always start in the background model

% When in the background state, we stay there with high prob
% When in the motif state, we immediately return to the background state.
transprob{1} = [0.8 0.2;
		1.0 0.0];


% LEVEL 2
startprob{2} = 'leftstart'; % both submodels start in substate 1
transprob{2} = zeros(motif_length, 2, motif_length);
termprob{2} = zeros(2, motif_length);

% In the background model, we only use state 1.
transprob{2}(1,1,1) = 1; % self loop
termprob{2}(1,1) = 0.2; % prob transition to end state

% Motif model
transprob{2}(:,2,:) = mk_leftright_transmat(motif_length, 0);
termprob{2}(2,end) = 1.0; % last state immediately terminates


% OBS LEVEl

obsprob = zeros([Qsize Osize]);
if 0
  % uniform background model
  obsprob(1,1,:) = normalise(ones(Osize,1));
else
  % deterministic background model (easy to see!)
  m = find(chars=='t');
  obsprob(1,1,m) = 1.0;
end
if 1
  % initialise with true motif (cheating)
  for i=1:motif_length
    m = find(chars == motif_pattern(i));
    obsprob(2,i,m) = 1.0;
  end
end

Oargs = {'CPT', obsprob};

[bnet, Qnodes, Fnodes, Onode] = mk_hhmm('Qsizes', Qsize, 'Osize', Osize, 'discrete_obs', 1, ...
	       'Oargs', Oargs, 'Ops', Qnodes(1:2), ...
	       'startprob', startprob, 'transprob', transprob, 'termprob', termprob);


Tmax = 20;
usecell = 0;

for seqi=1:5
  evidence = sample_dbn(bnet, Tmax, usecell);
  chars(evidence(end,:))
  %T = size(evidence, 2)
  %pretty_print_hhmm_parse(evidence, Qnodes, Fnodes, Onode, chars);
end
