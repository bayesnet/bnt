function bnet_miss = gener_MCAR_net(bnet_orig, base_proba)
% function bnet_miss = gener_MCAR_net(bnet_orig, base_proba)
% 
%   bnet_orig : a bnet
%   base_proba :  a probability for value to be missing
%
%   bnet_miss : a bnet that could be used in gener_data_from_bnet_miss function
%               to generate incomplete MCAR dataset
%
% Francois.Olivier.C.H@gmail.com

%%%%%%%%%%%% INIT
%bnet_miss = gener_MCAR_net(bnet_orig, base_proba, upd, bnet_biss)
%if nargin<4, upd = 0; else upd = 1; end
%if nargin<3, manual = 0; end
if nargin<2, error('Not enougth arguments'); end

% création du réseau
dag = bnet_orig.dag;
N = size(dag,2);
ns = bnet_orig.node_sizes;

%if nargin<4,
  ns_miss = zeros(1,3*N);
  ns_miss(1:N) = ns;
  ns_miss(N+1:2*N) = 2*ones(1,N); % 1= node i-N present, 2= node i-N missing
  ns_miss(2*N+1:3*N) = ns+1; % 1:ns, absent

  dag_miss = zeros(3*N,3*N);
  dag_miss(1:N,1:N) = dag;
  dag_miss(N+1:2*N,N+1:2*N)=mk_rnd_dag(N,N-ceil(rand*N/2)); %dag_miss(N+1:2*N,N+1:2*N)=mk_rnd_dag(N,2);
  for i=1:N, dag_miss(i,2*N+i)=1; dag_miss(N+i,2*N+i)=1; end

  bnet_miss = mk_bnet(dag_miss, ns_miss);
  CPT = CPT_from_bnet(bnet_orig, 0);
  for i=1:N
    bnet_miss.CPD{i} = tabular_CPD (bnet_miss, i, CPT{i});
  end

  % CPD of nodes M
  for i=1:N,
     if find(bnet_miss.order==i)<find(bnet_miss.order==N+i),
        CPT_M=[];
        for j=1:ns_miss(i), for l=1:ns_miss(N+i), for k=1:ns_miss(2*N+i),
          CPT_M=[CPT_M (((j==k)&(l==1))|((k==ns_miss(2*N+i))&(l==2)))];
        end, end, end
     else 
        CPT_M=[];
        for k=1:ns_miss(2*N+i), for l=1:ns_miss(N+i), for j=1:ns_miss(i), 
          CPT_M=[CPT_M (((j==k)&(l==1))|((k==ns_miss(2*N+i))&(l==2)))];
        end, end, end
     end
     bnet_miss.CPD{2*N+i} = tabular_CPD (bnet_miss, 2*N+i, CPT_M);
  end
%  else
%    ns_miss = bnet_miss.node_sizes;
%    dag_miss = bnet_miss.dag;
%  end

%%%%%%%%%%%% Base probability of missing value
p = base_proba;
 for i=1:N
  fam = find(dag_miss(:,N+i)==1)';
  semisize = prod(ns_miss(fam)); % as node N+i is binary to say i is present or missing
  CPT = zeros(1,2*semisize);
  CPT(1:semisize) = 1-p;
  CPT(semisize+1:2*semisize) = p;
  bnet_miss.CPD{N+i} = tabular_CPD (bnet_miss, N+i, CPT);
 end

%  if manual, % manual generation
%  
%  if upd==1;
%  
%   fprintf('Base probability of a data to be missing is %1.4f',base_proba);
%  
%   b=-1;
%   while ~(b==0 | b==1), b = input('Would you like to change a probability of a node to be missing (1 for yes, 0 for no) ?   '); end
%  
%   %%%%%%%%%%%% Update CPT with MCAR process
%   while b
%    fprintf('Nodes are from 1 to %d. ',N);
%    i=0;
%    while i<1 | i>N | round(i)~=i, i = input('Which node ?   '); end
%    fam = find(dag(:,i)==1)'; 
%    fam_miss = find(dag_miss(:,N+i)==1)';
%    cas = -ones(1, length(fam_miss)+1);
%    familly = [fam_miss, N+i];
%        p=-1;
%        while p<0 | p>1, p = input(' - A priori probability for this node to be missing ?   ');end
%        semisize = prod(ns_miss(fam_miss));
%        CPT = zeros(1,2*semisize);
%        CPT(1:semisize) = 1-p;
%        CPT(semisize+1:2*semisize) = p;
%        bnet_miss.CPD{N+i} = tabular_CPD (bnet_miss, N+i, CPT);
%    b=-1;
%    while ~(b==0 | b==1), b = input('Would you like to change a probability of a node to be missing (1 for yes, 0 for no) ?   ');end
%   end
%  end
%  
%  
%  else % automatic generation

%%%%%%%%%%%% creating BETAs
%      %% To use multiple of 5 percent in probs
%      if N<=25,
%          BETA = gener_problist(base_proba, N);
%      else
%          nboucles = floor(N/25);
%          BETA = [];
%          for i=1:nboucles
%              BETA1 = gener_problist(base_proba, 25);
%              BETA = [BETA, BETA1];
%          end
%          nreste = rem(N,25);
%          BETA1 = gener_problist(base_proba, nreste);
%          BETA = [BETA, BETA1];
%      end
   BETA = gener_discrete_dist(N, base_proba);

   order=[];
   missdagtmp = bnet_miss.dag(N+1:2*N,N+1:2*N);
   unprocessed = 1:N;
   while ~isempty(unprocessed)
     npar=[];
     for i=1:N, npar(end+1)=length(parents(missdagtmp,i));end,
     [npar, ord] = sort(npar);
     while ~ismember(ord(1),unprocessed)
       ord=ord(2:end);
     end
     order = [order, ord(1)];
     missdagtmp(ord(1),:)=0;
     unprocessed = mysetdiff(unprocessed,ord(1));
   end

%%%%%%%%%%%% Update CPT with MCAR process
for i=1:length(BETA)
  fam_miss = find(dag_miss(:,N+order(i))==1)';
  p=BETA(i);
  semisize = prod(ns_miss(fam_miss));

  if isempty(fam_miss),
        CPT = zeros(1,2*semisize);
        CPT(1:semisize) = 1-p;
        CPT(semisize+1:2*semisize) = p;
        bnet_miss.CPD{N+order(i)} = tabular_CPD (bnet_miss, N+order(i), CPT);
  else

    %node = N+order(i)
    %for k=1:semisize
    %  XI(k) = eval_xi(bnet_miss, N+order(i), k);
    %  %XI(k+semisize)=1-XI(k);
    %end
    %MUi1 = zeros(1,semisize);
    %MUi1 = gener_mu(p, semisize, XI);

        MUi1k = gener_discrete_dist(semisize, p);
        CPT = zeros(1,2*semisize);
        CPT(1:semisize) = 1-MUi1k;
        CPT(semisize+1:2*semisize) = MUi1k;
        bnet_miss.CPD{N+order(i)} = tabular_CPD (bnet_miss, N+order(i), CPT);

  end
% end
end