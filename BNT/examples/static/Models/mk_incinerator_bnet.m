function bnet  = mk_incinerator_bnet(ns)
% MK_INCINERATOR_BNET The waste incinerator emissions example from Cowell et al p145
% function bnet  = mk_incinerator_bnet(ns)
% 
% If ns is omitted, we use the scalars and binary nodes and the original params.
% Otherwise, we use random params of the desired size.
%
% Lauritzen, "Propogation of Probabilities, Means and Variances in Mixed Graphical Association Models", 
% JASA 87(420): 1098--1108
% This example is reprinted on p145 of "Probabilistic Networks and Expert Systems",
% Cowell, Dawid, Lauritzen and Spiegelhalter, 1999, Springer. 
% For a picture, see http://www.cs.berkeley.edu/~murphyk/Bayes/usage.html#cg_model

% node numbers
F = 1; W = 2; E = 3; B = 4; C = 5; D = 6; Min = 7; Mout = 8; L = 9;
names = {'F', 'W', 'E', 'B', 'C', 'D', 'Min', 'Mout', 'L'};
n = 9;
dnodes = [F W B];
cnodes = mysetdiff(1:n, dnodes);

% node sizes - all cts nodes are scalar, all discrete nodes are binary
if nargin < 1
  ns = ones(1, n);
  ns(dnodes) = 2;
  rnd = 0;
else
  rnd = 1;
end
  
% topology (p 1099, fig 1)
dag = zeros(n);
dag(F,E)=1;
dag(W,[E Min D]) = 1;
dag(E,D)=1;
dag(B,[C D])=1;
dag(D,[L Mout])=1;
dag(Min,Mout)=1;

% params (p 1102)
bnet = mk_bnet(dag, ns, 'discrete', dnodes, 'names', names);

if rnd
  for i=dnodes(:)'
    bnet.CPD{i} = tabular_CPD(bnet, i);
  end
  for i=cnodes(:)'
    bnet.CPD{i} = gaussian_CPD(bnet, i);
  end
else
  bnet.CPD{B} = tabular_CPD(bnet, B, 'CPT', [0.85 0.15]); % 1=stable, 2=unstable
  bnet.CPD{F} = tabular_CPD(bnet, F, 'CPT', [0.95 0.05]); % 1=intact, 2=defect
  bnet.CPD{W} = tabular_CPD(bnet, W, 'CPT', [2/7 5/7]); % 1=industrial, 2=household
  bnet.CPD{E} = gaussian_CPD(bnet, E, 'mean', [-3.9 -0.4 -3.2 -0.5], ...
			     'cov', [0.00002 0.0001 0.00002 0.0001]);
  bnet.CPD{D} = gaussian_CPD(bnet, D, 'mean', [6.5 6.0 7.5 7.0], ...
			     'cov', [0.03 0.04 0.1 0.1], 'weights', [1 1 1 1]);
  bnet.CPD{C} = gaussian_CPD(bnet, C, 'mean', [-2 -1], 'cov', [0.1 0.3]);
  bnet.CPD{L} = gaussian_CPD(bnet, L, 'mean', 3, 'cov', 0.25, 'weights', -0.5);
  bnet.CPD{Min} = gaussian_CPD(bnet, Min, 'mean', [0.5 -0.5], 'cov', [0.01 0.005]);
  bnet.CPD{Mout} = gaussian_CPD(bnet, Mout, 'mean', 0, 'cov', 0.002, 'weights', [1 1]);
end
