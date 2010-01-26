% Make the HHMM in Figure 1 of the NIPS'01 paper

Qsize = [2 3 2];
Qnodes = 1:3;
D = 3;
transprob = cell(1,D);
termprob = cell(1,D);
startprob = cell(1,D);
clear A;

% transprob{d}(i,k,j), transprob{1}(i,j)
% termprob{d}(k,j), termprob{1}(1,j)
% startprob{d}(k,j), startprob{1}(1,j)


% LEVEL 1

%       1 2 e
A{1} = [0 0 1;
	0 0 1];
[transprob{1}, termprob{1}] = remove_hhmm_end_state(A{1});
startprob{1} = [0.5 0.5];

% LEVEL 2
A{2} = zeros(Qsize(2), Qsize(1), Qsize(2)+1);

%              1 2 3 e
A{2}(:,1,:) = [0 1 0 0  % Q1=1 => model below state 0
	       0 0 1 0
	       0 0 0 1];

%              1 2 3 e
A{2}(:,2,:) = [0 1 0 0 % Q1=2 => model below state 1
	       0 0 1 0
	       0 0 0 1];

[transprob{2}, termprob{2}] = remove_hhmm_end_state(A{2});	       

% always enter level 2 in state 1
startprob{2} = [1 0 0
		1 0 0];

% LEVEL 3

A{3} = zeros([Qsize(3) Qsize(2) Qsize(3)+1]);
endstate = Qsize(3)+1;
%    Qt-1(3) Qt(2) Qt(3)
%                        1   2   e
A{3}(1,      1,    endstate) = 1.0; % Q2=1 => model below state 2/5
A{3}(:,      2,    :) = [0.0 1.0 0.0 % Q2=2 => model below state 3/6
            	         0.5 0.0 0.5]; 
A{3}(1,      3,    endstate) = 1.0; % Q2=3 => model below state 4/7

[transprob{3}, termprob{3}] = remove_hhmm_end_state(A{3});	       

startprob{3} = 'leftstart';



% OBS LEVEl

chars = ['a', 'b', 'c', 'd', 'x', 'y'];
Osize = length(chars);

obsprob = zeros([Qsize Osize]);
%       1 2 3 O
obsprob(1,1,1,find(chars == 'a')) =  1.0;

obsprob(1,2,1,find(chars == 'x')) =  1.0;
obsprob(1,2,2,find(chars == 'y')) =  1.0;

obsprob(1,3,1,find(chars == 'b')) =  1.0;

obsprob(2,1,1,find(chars == 'c')) =  1.0;

obsprob(2,2,1,find(chars == 'x')) =  1.0;
obsprob(2,2,2,find(chars == 'y')) =  1.0;

obsprob(2,3,1,find(chars == 'd')) =  1.0;

Oargs = {'CPT', obsprob};

bnet = mk_hhmm('Qsizes', Qsize, 'Osize', Osize, 'discrete_obs', 1, ...
	       'Oargs', Oargs, 'Ops', Qnodes(1:3), ...
	       'startprob', startprob, 'transprob', transprob, 'termprob', termprob);


Q1 = 1; Q2 = 2; Q3 = 3; F3 = 4; F2 = 5; Onode = 6;
Qnodes = [Q1 Q2 Q3]; Fnodes = [F2 F3];

for seqi=1:3
  evidence = sample_dbn(bnet, 'stop_test', 'is_F2_true_D3');      
  ev = cell2num(evidence);
  chars(ev(end,:))
  %T = size(evidence, 2)
  %pretty_print_hhmm_parse(evidence, Qnodes, Fnodes, Onode, chars);
end
