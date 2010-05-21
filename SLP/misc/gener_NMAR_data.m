function [data, comp_data, bnet_miss] = gener_NMAR_dataset(bnet_orig, m, bnet_miss, upd)
% [NMAR_data] = gener_NMAR_data(bnet_miss, length_of_dataset)
% 
% this function takes in input a bnet that can be used to generate NMAR data.
% this bnet bnet_miss can be creating by the function gener_NMAR_bnet.
%
% - comp_data (array) is a dataset that was generate by the bnet_orig that you enter in gener_NMAR_bnet
% - NMAR_data (cell array) is the dataset compdata that was emptyied by the NMAR process encodes in bnet_miss
%
% optional :
% - bnet_miss : an old bnet_miss built by this function
% - upd==1 if you want to update the bnet_miss
%
%  [data, comp_data, bnet_miss] = gener_NMAR_dataset(bnet_orig, m, bnet_miss, upd);
%
% version 0.5 : june 8th 2005, olivier.francois@insa-rouen.fr
%
% TO DO : 
%   - allow the combinaison of a missing state of one variable and another state of another variable to have influence
%   - allow the introduction of new nodes and specify which nodes it influence and which nodes has influence on it (it will also satisfy the first task then)
%


% INIT
N = size(bnet_orig.dag,2);
if mod(N,2)~=0, error('The number of nodes must be even'); end
if nargin<4, upd =0; end

% fisrt rules
if nargin<3,
 l1=[]; l2=[]; lp=[]; b=-1;
 while ~(b==0 | b==1), b = input('Would you like to make a node missing when another one is missing (1 for yes, 0 for no) ?   '); end
else
 b=-1;
 if upd, while ~(b==0 | b==1), b = input('Would you like to add rules (1 for yes, 0 for no) ?   '); end
   l1=bnet_miss.list{1};
   l2=bnet_miss.list{2};
   lp=bnet_miss.list{3}; 
 else
   l1=[]; l2=[]; lp=[]; b=-1;
 end
end
 
while b==1,
  n = input('The firts node ?   ');
  s = input('The node that have to be missing when this one is missing ?   ');
  p = input('The probability of the second node to be missing ?   ');  
  l1 = [l1, n]; l2 = [l2, s]; lp=[lp, p]; 
  b=-1;
  while ~(b==0 | b==1), b = input('Another one (1 for yes, 0 for no) ?   '); end
end
bb = length(l1);


if nargin>=3,
  bnet_miss = gener_NMAR_bnet(upd, bnet_orig, bnet_miss);
else
  bnet_miss = gener_NMAR_bnet(1, bnet_orig);
  bnet_miss.list={l1; l2; lp};
%%%%%%%%%%% SAVING FILE
 ss = 1;
 if nargin == 2 | upd==1, ss = input('Would you like to save the bnet of the NMAR process you have made (1 for yes) ?   '); end
 if ss == 1,
  ddd = datestr(now);
  ddd([12 15 18])='-' ;
  fnout=['NMAR-bnet-' ddd '.mat'];
  eval(['save ' fnout ' bnet_miss']); 
  fprintf(' The bnet for NMAR process was saved as : %s\n',fnout);
 end
end %if nargin

% Generation of a complete dataset
if N>9 & m>2000, disp('  ! It could take a long time...'); end
data = cell(N,m);
for l = 1:m, data(:,l) = sample_bnet(bnet_orig); end
disp('Complete data have been creating.');

% Generation of a NMAR dataset
miss_array = cell(2*N,m);
vide = cell(1,N); l= 1;
while l <= m, 
  ev(1:N) = data(:,l); ev(N+1:2*N) = vide;
  miss_array(:,l) = sample_bnet(bnet_miss, 'evidence', ev); 
  % apply simple rule of missingness
  ev2 = cell2mat(miss_array(N+1:2*N, l));
  if bb,
   missl1 = myintersect(find(ev2==2), l1);
   if ~isempty(missl1),
    for i=1:length(l1),
     if ev2(l1(i))==2, if rand<lp(i), ev2(l2(i))=2; miss_array{N+l2(i),l}=[2]; end, end
  end, end, end
  % verification that we have not a completly missing sample
  ev2 = 3-ev2; 
  if prod(ev2)==1, fprintf(' - %d, one completly missing sample removed', l); else l=l+1; end 
  if mod(l,100)==0, fprintf('\n - %d',l); end
end
fprintf('\n');
data = bnt_to_mat(data); comp_data = data;
miss_array = bnt_to_mat(miss_array(N+1:2*N, :));
miss_array = 2-miss_array;
data = data.*miss_array;
data = mat_to_bnt(data, 0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bnet_miss = gener_NMAR_bnet(upd, bnet_orig, bnet_miss)

%%%%%%%%%%%% INIT
dag = bnet_orig.dag;
N = size(dag,2);
ns = bnet_orig.node_sizes;
NN = 0;
ns_miss = zeros(1,2*N+NN);
ns_miss(1:N) = ns;
ns_miss(N+1:2*N) = 2*ones(1,N); % 1= node i-N present, 2= node i-N missing
  dagm = zeros(2*N,2*N);
  dagm(1:N,1:N) = dag;
  dagm(N+1:2*N, N+1:2*N) = dag;
  dagm(1:N, N+1:2*N) = dag;
  for i=1:N, dagm(i, i+N)=1; end

%%%%%%%%%%%% NEW NODES
%  b=-1;
%  while ~(b==0 | b==1), b = input('Would you like to add new nodes ?   '); end
%  
%  if b==1,
%  NN = -1;
%  while (NN<0 | round(NN)~=NN), NN = input('How Many ?   '); end
%  if NN==1, fprintf('New node will be called %d\n',(2*N+1));
%  else fprintf('New node will be called %d and following numbers\n',(2*N+1)); end
%  
%  disp(' !!! Make sur that the dependence you will create will not create cycle in the Bnet used to create NMAR data !!! here is the current DAG');
%  dagm2 = zeros(2*N+NN, 2*N+NN);
  dagm2(1:2*N, 1:2*N) = dagm; 
%  draw_graph(dagm2); drawnow;
%  clear dagm
%  
%  L1={};L2={};
%  for i=1:NN
%    fprintf('For the node %d, ',(2*N+i));
%    L1{i} = input('FROM which nodes will it have influence (sample [3 2]) ?   ');
%    L2{i} = input('ON which nodes will it have influence (sample [1 4 3]) ?   '); 
%    ns_miss(2*N+i) = input('What is its size ?   '); 
%  end
%  
%  end

%%%%%%%%%%%% BNET CREATION
if nargin==2,
  dag_miss = dagm2;
  for i=1:NN, dag_miss(L1{i},2*N+i)=1; dag_miss(2*N+i, L2{i})=1; end

  bnet_miss = mk_bnet(dag_miss, ns_miss);
  CPT = CPT_from_bnet(bnet_orig);
  for i=1:N
    bnet_miss.CPD{i} = tabular_CPD (bnet_miss, i, CPT{i}); % error with new nodes
  end
elseif nargin==3,
  ns_miss = bnet_miss.node_sizes;
  dag_miss = bnet_miss.dag;
end
clear dagm2
order = 1:(2*N+NN);

%%%%%%%%%%%% Base probability of missing value
if nargin==2, b=1; else b=0; end
if b==1,
disp('Probability MUST be between 0 and 1.');

if nargin==2, 
 p=-1;
 while p<0 | p>1, p = input('Base probability of a value to be missing ?   '); end
 for i=1:N
  fam = find(dag_miss(:,N+i)==1)';
  semisize = prod(ns_miss(fam)); % as node N+i is binary to say i is present or missing
  CPT = zeros(1,2*semisize);
  CPT(1:semisize) = 1-p;
  CPT(semisize+1:2*semisize) = p;
  bnet_miss.CPD{N+i} = tabular_CPD (bnet_miss, N+i, CPT);
 end
end
end

b=-1;
while ~(b==0 | b==1), b = input('Would you like to change a probability of a node to be missing (1 for yes, 0 for no) ?   '); end
if b, disp(' BE CAREFULL !! New rules can overwrite old ones partialy or fully !! So the order of entries is important'); end

%%%%%%%%%%%% Update CPT with NMAR process
while b
  fprintf('Nodes are from 1 to %d. ',N);
  i=0;
  while i<1 | i>N | round(i)~=i, i = input('Which node ?   '); end
  fam = find(dag(:,i)==1)'; 
  fam_miss = find(dag_miss(:,N+i)==1)';
  cas = -ones(1, length(fam_miss)+1);
  familly = [fam_miss, N+i];
  fprintf('States are from 1 to %d (-1 for any states, -2 to cancel). For which state of the variable %d ?', ns(i), i);
  state=-3;
  while state<-2 | state>ns(i) | round(state)~=state | state==0, state  = input('   ');end
  if isempty(fam),
    if state==-1,
      p=-1;
      while p<0 | p>1, p = input(' - A priori probability for this node to be missing ?   ');end
      semisize = prod(ns_miss(fam_miss));
      CPT = zeros(1,2*semisize);
      CPT(1:semisize) = 1-p;
      CPT(semisize+1:2*semisize) = p;
      bnet_miss.CPD{N+i} = tabular_CPD (bnet_miss, N+i, CPT);
    elseif state~=-2
      cas = state;
      CPT = CPT_from_bnet(bnet_miss);
      CPT = CPT{N+i};
      p=-1;
      while p<0 | p>1, p = input(' - A priori probability for this node to be missing in this state ?   ');end
      ind = subv2ind(ns_miss(familly),[cas, 1]);
      CPT(ind)=1-p;
      ind = subv2ind(ns_miss(familly),[cas, 2]);
      CPT(ind)=p;
      bnet_miss.CPD{N+i} = tabular_CPD (bnet_miss, N+i, CPT);      
    end
  else 
   if state>-2,
    siz=length(cas);
    place = find(fam_miss==i);
    cas(place) = state;
    
    for k = fam, 
        state=-3;
        fprintf(' - For the parent named %d, states are from 1 to %d (-1 for any states of this parent). ',k, ns(k));
        while state<0 | state>ns(k) | round(state)~=state, state  = input('Which state ?   ');end
        %if state==0,
        %  place = find(fam_miss==(fam_miss(k)+N));
        %  cas(place) = 2;                   % Missing
        %elseif state==-2,                   % a changer ???
        %  disp(' This case is buggy, taking missing state instand to minimise influence.');
        %  place = find(fam_miss==(fam_miss(k)+N));
        %  cas(place) = 2;  
        %elseif state~=0 & state~=-2, 
          place = find(fam_miss==(fam_miss(k)));
          cas(place) = state;               % Present
          if state~=-1; place = find(fam_miss==(fam_miss(k)+N)); cas(place) = 1; end 
        %end
      end
    
    p=-1;
    while p<0 | p>1, p = input('Probability in this case of the value to be missing ?   ');end
    CPT = CPT_from_bnet(bnet_miss);
    CPT = CPT{N+i};
    
    cas(end) = 1; % i is present
    subcas_names = find(cas==-1);
    if isempty(subcas_names),
      ind = subv2ind(ns_miss(familly),cas);
      CPT(ind) = 1-p;
    else
      subcas = ones(1, length(subcas_names));
      continu = 1;
      while continu
        cas(subcas_names) = subcas;
        ind = subv2ind(ns_miss(familly),cas);
        CPT(ind) = 1-p;
        [subcas, continu] = next_case(subcas, ns_miss(familly(subcas_names)));
      end
    end
      
    cas(end)=2; % i is missing
    if isempty(subcas_names),
      ind = subv2ind(ns_miss(familly),cas);
      CPT(ind) = p;
    else
      subcas = ones(1, length(subcas_names));
      continu = 1;
      while continu
        cas(subcas_names) = subcas;
        ind = subv2ind(ns_miss(familly),cas);
        CPT(ind) = p;
        [subcas, continu] = next_case(subcas, ns_miss(familly(subcas_names)));
      end
    end
    mass = sum(CPT, length(size(CPT)));
    while length(size(mass))>2, mass = prod(mass, length(size(mass))); end
    mass = prod(prod(mass));
    if mass~=1, disp('not a proba...'); end
    bnet_miss.CPD{N+i} = tabular_CPD (bnet_miss, N+i, CPT);
   end %if state~=-2 for the node
  end %if isempty(fam),
  b=-1;
  while ~(b==0 | b==1), b = input('Would you like to change a probability of a node to be missing (1 for yes, 0 for no) ?   ');end
end

