function sample = sample_bnet(bnet, varargin)
% SAMPLE_BNET Generate a random sample from a Bayes net.
% SAMPLE = SAMPLE_BNET(BNET, ...)
%
% sample{i} contains the value of the i'th node.
% i.e., the result is an Nx1 cell array.
% Nodes are sampled in the order given by bnet.order.
%
% Optional arguments:
%
% evidence - initial evidence; if evidence{i} is non-empty, node i won't be sampled.

% set defauly params
n = length(bnet.dag);
sample = cell(n,1);

% get optional params
args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'evidence',    sample = args{i+1}(:);
   otherwise, error(['unrecognized argument ' args{i}])
  end
end

for j=bnet.order(:)'
  if isempty(sample{j})
    %ps = parents(bnet.dag, j);
    ps = bnet.parents{j};
    e = bnet.equiv_class(j);
    sample{j} = sample_node(bnet.CPD{e}, sample(ps));
  end
end
