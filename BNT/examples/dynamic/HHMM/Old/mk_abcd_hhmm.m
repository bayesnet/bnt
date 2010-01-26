% Make the HHMM in Figure 1 of the NIPS'01 paper

Qsize = [2 3 2];
D = 3;

% transprob{d}(i,k,j), transprob{1}(i,j)
% termprob{d}(k,j), termprob{1}(1,j)
% startprob{d}(k,j), startprob{1}(1,j)
% obsprob(k, o) for discrete outputs
    
% LEVEL 1
%       1 2 e
A{1} = [0 0 1;
	0 0 1];
[transprob{1}, termprob{1}] = remove_hhmm_end_state(A{1});
startprob{1} = [0.5 0.5];
Q1args = {'startprob', startprob{1}, 'transprob', transprob{1}};

% LEVEL 2
A{2} = zeros(Qsize(2), Qsize(1), Qsize(2)+1);

%              1 2 3 e
A{2}(:,1,:) = [0 1 0 0
	       0 0 1 0
	       0 0 0 1];

%              1 2 3 e
A{2}(:,2,:) = [0 1 0 0
	       0 0 1 0
	       0 0 0 1];

[transprob{2}, termprob{2}] = remove_hhmm_end_state(A{2});	       

% always enter level 2 in state 1
startprob{2} = [1 0 0
		1 0 0];

Q2args = {'startprob', startprob{2}, 'transprob', transprob{2}};
F2args = {'CPT', termprob{2}};


% LEVEL 3

A{3} = zeros([Qsize(3) Qsize(1:2) Qsize(3)+1]);
endstate = Qsize(3)+1;
%    Qt-1(3) Qt(1) Qt(2) Qt(3)
%                               1   2   e
A{3}(1,      1,    1,    endstate) = 1.0;
A{3}(:,      1,    2,    :) = [0.0 1.0 0.0
            	               0.5 0.0 0.5];
A{3}(1,      1,    3,    endstate) = 1.0;

A{3}(1,      2,    1,    endstate) = 1.0;
A{3}(:,      2,    2,    :) = [0.0 1.0 0.0
            	               0.5 0.0 0.5];
A{3}(1,      2,    3,    endstate) = 1.0;

A{3} = reshape(A{3}, [Qsize(3) prod(Qsize(1:2)) Qsize(3)+1]);
[transprob{3}, termprob{3}] = remove_hhmm_end_state(A{3});	       

% define the vertical entry points to level 3
startprob{3} = zeros(Qsize);
%            Q1 Q2 Q3
startprob{3}(1, 1, 1) = 1.0;
startprob{3}(1, 2, 1) = 1.0;
startprob{3}(1, 3, 1) = 1.0;

startprob{3}(2, 1, 1) = 1.0;
startprob{3}(2, 2, 1) = 1.0;
startprob{3}(2, 3, 1) = 1.0;

startprob{3} = reshape(startprob{3}, prod(Qsize(1:2)), Qsize(3));

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

obsprob = reshape(obsprob, prod(Qsize), Osize);

[intra, inter, Qnodes, Fnodes, Onode] = mk_hhmm_topo(D);

hhmm.Qnodes = Qnodes;
hhmm.Fnodes = Fnodes;
hhmm.Onode = Onode;
hhmm.D = D;
hhmm.Qsize = Qsize;
hhmm.Osize = Osize;
hhmm.startprob = startprob;
hhmm.transprob = transprob;
hhmm.termprob = termprob;
hhmm.obsprob = obsprob;
hhmm.A = A;


