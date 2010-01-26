function bnet = mk_chmm(N, Q, Y, discrete_obs, coupled, CPD)
% MK_CHMM Make a coupled Hidden Markov Model
%
% There are N hidden nodes, each connected to itself and its two nearest neighbors in the next
% slice (apart from the edges, where there is 1 nearest neighbor).
%
% Example: If N = 3, the hidden backbone is as follows, where all arrows point to the righ+t
%
% X1--X2
%   \/ 
%   /\
% X2--X2
%   \/ 
%   /\
% X3--X3
%
% Each hidden node has a "private" observed child (not shown).
%
% BNET = MK_CHMM(N, Q, Y)
% Each hidden node is discrete and has Q values.
% Each observed node is a Gaussian vector of length Y.
%
% BNET = MK_CHMM(N, Q, Y, DISCRETE_OBS)
% If discrete_obs = 1, the observations are discrete (values in {1, .., Y}).
%
% BNET = MK_CHMM(N, Q, Y, DISCRETE_OBS, COUPLED)
% If coupled = 0, the chains are not coupled, i.e., we make N parallel HMMs.
%
% BNET = MK_CHMM(N, Q, Y, DISCRETE_OBS, COUPLED, CPDs)
% means use the specified CPD structures instead of creating random params.
%  CPD{i}.CPT, i=1:N specifies the prior
%  CPD{i}.CPT, i=2N+1:3N specifies the transition model
%  CPD{i}.mean, CPD{i}.cov, i=N+1:2N specifies the observation model if Gaussian
%  CPD{i}.CPT, i=N+1:2N if discrete


if nargin < 2, Q = 2; end
if nargin < 3, Y = 1; end
if nargin < 4, discrete_obs = 0; end
if nargin < 5, coupled = 1; end
if nargin < 6, rnd = 1; else rnd = 0; end
  
ss = N*2;
hnodes = 1:N;
onodes = (1:N)+N;

intra = zeros(ss);
for i=1:N
  intra(hnodes(i), onodes(i))=1;
end

inter = zeros(ss);
if coupled
  for i=1:N
    inter(i, max(i-1,1):min(i+1,N))=1;
  end
else
  inter(1:N, 1:N) = eye(N);
end  

ns = [Q*ones(1,N) Y*ones(1,N)]; 

eclass1 = [hnodes onodes];
eclass2 = [hnodes+ss onodes];
if discrete_obs
  dnodes = 1:ss;
else
  dnodes = hnodes;
end
bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'eclass1', eclass1, 'eclass2', eclass2, ...
	      'observed', onodes);

if rnd
  for i=hnodes(:)'
    bnet.CPD{i} = tabular_CPD(bnet, i);
  end
  for i=onodes(:)'
    if discrete_obs
      bnet.CPD{i} = tabular_CPD(bnet, i);
    else
      bnet.CPD{i} = gaussian_CPD(bnet, i);
    end
  end
  for i=hnodes(:)'+ss
    bnet.CPD{i} = tabular_CPD(bnet, i);
  end
else
  for i=hnodes(:)'
    bnet.CPD{i} = tabular_CPD(bnet, i, CPD{i}.CPT);
  end
  for i=onodes(:)'
    if discrete_obs
      bnet.CPD{i} = tabular_CPD(bnet, i, CPD{i}.CPT);
    else
      bnet.CPD{i} = gaussian_CPD(bnet, i, CPD{i}.mean, CPD{i}.cov);
    end
  end
  for i=hnodes(:)'+ss
    bnet.CPD{i} = tabular_CPD(bnet, i, CPD{i}.CPT);
  end
end


