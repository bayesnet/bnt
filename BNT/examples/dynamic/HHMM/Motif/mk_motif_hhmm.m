function [bnet, Qnodes, Fnodes, Onode] = mk_motif_hhmm(varargin)
% [bnet, Qnodes, Fnodes, Onode] = mk_motif_hhmm(...)
%
% Make the following HHMM
%
%    S2 <----------------------> S1
%    |                           |
%    |                           |
%   M1 -> M2 -> M3 -> end        B1 -> end
%
% where Mi represents the i'th letter in the motif
% and B is the background state.
% Si chooses between running the motif or the background.
% The Si and B states have self loops (not shown).
%
% The transition params are defined to respect the above topology.
% The background is uniform; each motif state has a random obs. distribution.
%
% Optional params:
% motif_length  - required, unless we specify motif_pattern
% motif_pattern - if specified, we make the motif submodel deterministically
%                  emit this pattern
% background    - if specified, we make the background submodel
%                  deterministically emit this (makes the motif easier to see!)


args = varargin;
nargs = length(args);

% extract pattern, if any
motif_pattern = [];
for i=1:2:nargs
  switch args{i},
   case 'motif_pattern', motif_pattern = args{i+1}; 
  end
end

% set defaults
motif_length = length(motif_pattern);
background_char = [];

% get params
for i=1:2:nargs
  switch args{i},
   case 'motif_length', motif_length = args{i+1}; 
   case 'background', background_char = args{i+1};
  end
end


chars = ['a', 'c', 'g', 't'];
Osize = length(chars);

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
transprob{2}(:,2,:) = mk_leftright_transmat(motif_length, 0); % no self loops
termprob{2}(2,end) = 1.0; % last state immediately terminates


% OBS LEVEl

obsprob = zeros([Qsize Osize]);
if isempty(background_char)
  % uniform background model
  %obsprob(1,1,:) = normalise(ones(Osize,1));
  obsprob(1,1,:) = normalise(rand(Osize,1));
else
  % deterministic background model (easy to see!)
  m = find(chars==background_char);
  obsprob(1,1,m) = 1.0;
end

if ~isempty(motif_pattern)
  % initialise with true motif (cheating)
  for i=1:motif_length
    m = find(chars == motif_pattern(i));
    obsprob(2,i,m) = 1.0;
  end
else
  obsprob(2,:,:) = mk_stochastic(rand(motif_length, Osize));
end

if 0
  Oargs = {'CPT', obsprob};
else
  % We use a minent prior for the emission distribution for the states in the motif model
  % (but not the background model). This encourages nearly deterministic distributions.
  % We create an index matrix  (where M = motif length)
  %  [2 1
  %   2 2
  %   ...
  %   2 M]
  % and then convert this to a list of integers, which
  % specifies when to use the minent prior (Q1=2 specifies motif model).
  M = motif_length;
  ndx = [2*ones(M,1) (1:M)'];
  pcases = subv2ind([2 motif_length], ndx);
  Oargs = {'CPT', obsprob, 'prior_type', 'entropic', 'entropic_pcases', pcases};
end



[bnet, Qnodes, Fnodes, Onode] = mk_hhmm('Qsizes', Qsize, 'Osize', Osize, 'discrete_obs', 1, ...
	       'Oargs', Oargs, 'Ops', Qnodes(1:2), ...
	       'startprob', startprob, 'transprob', transprob, 'termprob', termprob);

