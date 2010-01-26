function bnet = mk_fhmm(N, Q, Y, discrete_obs)
% MK_FHMM Make a factorial Hidden Markov Model
%
% There are N independent parallel hidden chains, each connected to the output
%
% e.g., N = 2 (vertical/diagonal edges point down)
%
% A1--->A2
% | B1--|->B2
% | /   |/
% Y1    Y2
%
% [bnet, onode] = mk_chmm(n, q, y, discrete_obs)
%
% Each hidden node is discrete and has Q values.
% If discrete_obs = 1, each observed node is discrete and has values 1..Y.
% If discrete_obs = 0, each observed node is a Gaussian vector of length Y.

if nargin < 2, Q = 2; end
if nargin < 3, Y = 2; end
if nargin < 4, discrete_obs = 1; end

ss = N+1;
hnodes = 1:N;
onode = N+1;

intra = zeros(ss);
intra(hnodes, onode) = 1;

inter = eye(ss);
inter(onode,onode) = 0;

ns = [Q*ones(1,N) Y];

eclass1 = [hnodes onode];
eclass2 = [hnodes+ss onode];
if discrete_obs
  dnodes = 1:ss;
else
  dnodes = hnodes;
end
bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'eclass1', eclass1, 'eclass2', eclass2, ...
	      'observed', onode);

for i=hnodes(:)'
  bnet.CPD{i} = tabular_CPD(bnet, i);
end
i = onode;
if discrete_obs
  bnet.CPD{i} = tabular_CPD(bnet, i);
else
  bnet.CPD{i} = gaussian_CPD(bnet, i);
end
for i=hnodes(:)'+ss
  bnet.CPD{i} = tabular_CPD(bnet, i);
end


