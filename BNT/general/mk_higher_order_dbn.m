function bnet = mk_higher_order_dbn(intra, inter, node_sizes, varargin)
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

% As this method is used to generate a higher order Markov Model
% also connect from time slice t - i -> t with i > 1 has to be 
% taken into account.

%inter should be a three dimensional array where inter(:,:,i)
%describes the connections from time-slice t - i to t.  
[rows,columns,order] = size(inter);
assert(rows    == n);
assert(columns == n);
dag = zeros((order + 1)*n);

i = 0;
while i <= order
    j = i;
    while j <= order
        if j == i
            dag(1 + i*n:(i+1)*n,1+i*n:(i+1)*n) = intra;
        else
            dag(1+i*n:(i+1)*n,1+j*n:(j+1)*n) = inter(:,:,j - i);
        end
        j = j + 1;
    end;
    i = i + 1;
end;

bnet.dag = dag;
bnet.names = {};

directed = 1;
if ~acyclic(dag,directed)
  error('graph must be acyclic')
end

% Calculation of the equivalence classes
bnet.eclass1 = 1:n;
bnet.eclass = zeros(order + 1,ss);
bnet.eclass(1,:) = 1:n;
for i = 1:order
    bnet.eclass(i+1,:) = bnet.eclass(i,:);
    for j = 1:ss 
        if(isequal(parents(dag,(i-1)*n+j)+ss,parents(dag,(i*n + j))))
	   %fprintf('%d has isomorphic parents, eclass %d \n',j,bnet.eclass(i,j))
        else
	   bnet.eclass(i + 1,j) = max(bnet.eclass(i+1,:))+1;
	   %fprintf('%d has non isomorphic parents, eclass %d \n',j,bnet.eclass(i,j))  
	end;
    end;
end;
bnet.eclass1 = 1:n;

% To be compatible with whe rest of the code 
bnet.eclass2 = bnet.eclass(2,:);

dnodes = 1:n;
bnet.observed = [];

if nargin >= 4
  args = varargin;
  nargs = length(args);
  if ~isstr(args{1})
    if nargs >= 1 dnodes = args{1}; end
    if nargs >= 2 bnet.eclass1 = args{2}; bnet.eclass(1,:) = args{2}; end
    if nargs >= 3 bnet.eclass2 = args{3}; bnet.eclass(2,:) = args{2}; end
    if nargs >= 4 bnet.intra1 = args{4}; end
  else
    for i=1:2:nargs
      switch args{i},
       case 'discrete', dnodes = args{i+1}; 
       case 'observed', bnet.observed = args{i+1}; 
       case 'eclass1',  bnet.eclass1 = args{i+1}; bnet.eclass(1,:) = args{i+1}; 
       case 'eclass2',  bnet.eclass2 = args{i+1}; bnet.eclass(2,:) = args{i+1};
       case 'eclass',   bnet.eclass = args{i+1};  
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
bnet.node_sizes = repmat(ns(:),1,order + 1);

cnodes = mysetdiff(1:n, dnodes);
bnet.dnodes_slice = dnodes;
bnet.cnodes_slice = cnodes;
bnet.dnodes = dnodes;
bnet.cnodes = cnodes;
% To adapt the function to higher order Markov models include dnodes for more 
% time slices
for i = 1:order
    bnet.dnodes = [bnet.dnodes dnodes+i*n];
    bnet.cnodes = [bnet.cnodes cnodes+i*n];
end

% Generieren einer Matrix, deren i-te Spalte die Aequivalenzklassen
% der i-ten Zeitscheibe enthaelt. 
bnet.equiv_class = [bnet.eclass(1,:)]';
for i = 2:(order + 1)
    bnet.equiv_class = [bnet.equiv_class   bnet.eclass(i,:)'];
end

bnet.CPD = cell(1,max(bnet.equiv_class(:)));

ss = n;
onodes = bnet.observed;
hnodes = mysetdiff(1:ss, onodes);
bnet.hidden_bitv = zeros(1,(order + 1)*ss);
for i = 0:order
    bnet.hidden_bitv(hnodes +i*ss) = 1;
end;

bnet.parents = cell(1, (order + 1)*ss);
for i=1:(order + 1)*ss
  bnet.parents{i} = parents(bnet.dag, i);
end

bnet.auto_regressive = zeros(1,ss);
% ar(i)=1 means (observed) node i depends on i in the  previous slice
for o=bnet.observed(:)'
  if any(bnet.parents{o+ss} <= ss)
    bnet.auto_regressive(o) = 1;
  end
end















