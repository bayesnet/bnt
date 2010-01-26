function bnet = mk_dbn(intra, inter, node_sizes, varargin)
% MK_DBN Make a Dynamic Bayesian Network.
%
% BNET = MK_DBN(INTRA, INTER, NODE_SIZES, ...) makes a DBN with arcs
% from i in slice t to j in slice t iff intra(i,j) = 1, and 
% from i in slice t to j in slice t+1 iff inter(i,j) = 1,
% for i,j in {1, 2, ..., n}, where n = num. nodes per slice, and t >= 1.
% node_sizes(i) is the number of values node i can take on.
% The nodes are assumed to be in topological order. Use TOPOLOGICAL_SORT if necessary.
% See also mk_bnet.
%
% Optional arguments [default in brackets]
% 'discrete' - list of discrete nodes [1:n]
% 'observed' - the list of nodes which will definitely be observed in every slice of every case [ [] ]
% 'eclass1' - equiv class for slice 1 [1:n]
% 'eclass2' - equiv class for slice 2 [tie nodes with equivalent parents to slice 1]
%    equiv_class1(i) = j means node i in slice 1 gets its parameters from bnet.CPD{j},
%    i.e., nodes i and j have tied parameters.
% 'intra1' - topology of first slice, if different from others
% 'names' - a cell array of strings to be associated with nodes 1:n [{}]
%    This creates an associative array, so you write e.g.
%     'evidence(bnet.names{'bar'}) = 42' instead of  'evidence(2} = 42' 
%     assuming names = { 'foo', 'bar', ...}.
%    
% For backwards compatibility with BNT2, arguments can also be specified as follows
%   bnet = mk_dbn(intra, inter, node_sizes, dnodes, eclass1, eclass2, intra1)
%
% After calling this function, you must specify the parameters (conditional probability
% distributions) using bnet.CPD{i} = gaussian_CPD(...) or tabular_CPD(...) etc.


n = length(intra);
ss = n;
bnet.nnodes_per_slice = ss;
bnet.intra = intra;
bnet.inter = inter;
bnet.intra1 = intra;
dag = zeros(2*n);
dag(1:n,1:n) = bnet.intra1;
dag(1:n,(1:n)+n) = bnet.inter;
dag((1:n)+n,(1:n)+n) = bnet.intra;
bnet.dag = dag;
bnet.names = {};

directed = 1;
if ~acyclic(dag,directed)
  error('graph must be acyclic')
end


bnet.eclass1 = 1:n;
%bnet.eclass2 = (1:n)+n;
bnet.eclass2 = bnet.eclass1;
for i=1:ss
  if isequal(parents(dag, i+ss), parents(dag, i)+ss)
    %fprintf('%d has isomorphic parents, eclass %d\n', i, bnet.eclass2(i))
  else
    bnet.eclass2(i) = max(bnet.eclass2) + 1;
    %fprintf('%d has non isomorphic parents, eclass %d\n', i, bnet.eclass2(i))
  end
end

dnodes = 1:n;
bnet.observed = [];

if nargin >= 4
  args = varargin;
  nargs = length(args);
  if ~isstr(args{1})
    if nargs >= 1, dnodes = args{1}; end
    if nargs >= 2, bnet.eclass1 = args{2}; end
    if nargs >= 3, bnet.eclass2 = args{3}; end
    if nargs >= 4, bnet.intra1 = args{4}; end
  else
    for i=1:2:nargs
      switch args{i},
       case 'discrete', dnodes = args{i+1}; 
       case 'observed', bnet.observed = args{i+1}; 
       case 'eclass1',  bnet.eclass1 = args{i+1}; 
       case 'eclass2',  bnet.eclass2 = args{i+1}; 
       case 'intra1',  bnet.intra1 = args{i+1}; 
       %case 'ar_hmm',  bnet.ar_hmm = args{i+1};  % should check topology
       case 'names',  bnet.names = assocarray(args{i+1}, num2cell(1:n)); 
       otherwise,  
	error(['invalid argument name ' args{i}]);       
      end
    end
  end
end


bnet.observed = sort(bnet.observed); % for comparing sets
ns = node_sizes;
bnet.node_sizes_slice = ns(:)';
bnet.node_sizes = [ns(:) ns(:)];

cnodes = mysetdiff(1:n, dnodes);
bnet.dnodes_slice = dnodes;
bnet.cnodes_slice = cnodes;
bnet.dnodes = [dnodes dnodes+n];
bnet.cnodes = [cnodes cnodes+n];

bnet.equiv_class = [bnet.eclass1(:) bnet.eclass2(:)];
bnet.CPD = cell(1,max(bnet.equiv_class(:)));
eclass = bnet.equiv_class(:);
E = max(eclass);
bnet.rep_of_eclass = zeros(1,E);
for e=1:E
  mems = find(eclass==e);
  bnet.rep_of_eclass(e) = mems(1);
end

ss = n;
onodes = bnet.observed;
hnodes = mysetdiff(1:ss, onodes);
bnet.hidden_bitv = zeros(1,2*ss);
bnet.hidden_bitv(hnodes) = 1;
bnet.hidden_bitv(hnodes+ss) = 1;

bnet.parents = cell(1, 2*ss);
for i=1:ss
  bnet.parents{i} = parents(bnet.dag, i);
  bnet.parents{i+ss} = parents(bnet.dag, i+ss);
end

bnet.auto_regressive = zeros(1,ss);
% ar(i)=1 means (observed) node i depends on i in the  previous slice
for o=bnet.observed(:)'
  if any(bnet.parents{o+ss} <= ss)
    bnet.auto_regressive(o) = 1;
  end
end

