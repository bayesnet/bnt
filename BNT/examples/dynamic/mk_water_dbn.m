function bnet = mk_water_dbn(discrete_obs, obs_leaves)
% MK_WATER_DBN
% bnet = mk_water_dbn(discrete_obs, obs_leaves)
%
% If discrete_obs = 1 (default), the leaves are binary, else scalar Gaussians
% If obs_leaves = 1, all the leaves are observed, otherwise rnd nodes are observed
%
% This is a model of the biological processes of a water purification plant, developed
% by Finn V. Jensen, Uffe Kjærulff, Kristian G. Olesen, and Jan Pedersen.
% See http://www-nt.cs.berkeley.edu/home/nir/public_html/Repository/water.htm
% See also Boyen and Koller, "Tractable Inference for Complex Stochastic Processes", UAI98

if nargin < 1, discrete_obs = 1; end
if nargin < 1, obs_leaves = 1; end

ss = 12;
intra = zeros(ss);
intra(1,9) = 1;
intra(3,10) = 1;
intra(4,11) = 1;
intra(8,12) = 1;

inter = zeros(ss);
inter(1, [1 3]) = 1;
inter(2, [2 3 7]) = 1;
inter(3, [3 4 5]) = 1;
inter(4, [3 4 6]) = 1;
inter(5, [3 5 6]) = 1;
inter(6, [4 5 6]) = 1;
inter(7, [7 8]) = 1;
inter(8, [6 7 8]) = 1;

if obs_leaves
  onodes = 9:12; % leaves
else
  onodes = [1 5 9:12]; % throw in some other nodes
end
hnodes = 1:8;
if discrete_obs
  ns = 2*ones(1 ,ss);
  dnodes = 1:ss;
else
  ns = [2*ones(1,length(hnodes)) 1*ones(length(onodes))];
  dnodes = hnodes;
end

eclass1 = 1:12;
eclass2 = [13:20 9:12];
bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'eclass1', eclass1, 'eclass2', eclass2, ...
	      'observed', onodes);
if discrete_obs
  for i=1:max(eclass2)
    bnet.CPD{i} = tabular_CPD(bnet, i);
  end
else
  for i=hnodes(:)'
    bnet.CPD{i} = tabular_CPD(bnet, i);
  end
  for i=onodes(:)'
    bnet.CPD{i} = gaussian_CPD(bnet, i);
  end
  for i=hnodes(:)'+ss
    bnet.CPD{i} = tabular_CPD(bnet, i);
  end
end



