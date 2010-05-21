function bnet_miss = gener_MAR_net(bnet_orig, base_proba)
% function bnet_miss = gener_MAR_net(bnet_orig, base_proba)
% 
%   bnet_orig : a bnet
%   base_proba :  a probability for value to be missing
%
%   bnet_miss : a bnet that could be used in gener_data_from_bnet_miss function
%               to generate incomplete MAR dataset
%
% Francois.Olivier.C.H@gmail.com

%%%%%%%%%%%% INIT
if nargin<2, error('Not enougth arguments'); end

% création du réseau
dag = bnet_orig.dag;
N = size(dag,2);
ns = bnet_orig.node_sizes;

  ns_miss = zeros(1,3*N);
  ns_miss(1:N) = ns;
  ns_miss(N+1:2*N) = 2*ones(1,N); % 1= node i-N present, 2= node i-N missing
  ns_miss(2*N+1:3*N) = ns+1; % 1:ns, absent

  dag_miss = zeros(3*N,3*N);
  dag_miss(1:N,1:N) = dag;
%  dag_miss(2*N+1:3*N,N+1:2*N)=mk_rnd_dag(N,N-ceil(rand*N/2)); 
  lim = 1+(rand>.4)+(rand>.65)+(rand>.9);
  dag_miss(2*N+1:3*N,N+1:2*N)=mk_rnd_dag(N,lim); 

  for i=1:N, dag_miss(i,2*N+i)=1; dag_miss(N+i,2*N+i)=1; dag_miss(2*N+i,i)=0; dag_miss(2*N+i,N+i)=0; end

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

   BETA = gener_discrete_dist(N, base_proba);

   order=[];
   missdagtmp = bnet_miss.dag; %(N+1:2*N,N+1:2*N);
   unprocessed = 1:N;
   while ~isempty(unprocessed)
     npar=[];
     for i=N+1:2*N, npar(end+1)=length(parents(missdagtmp,i));end,  % to be verifie from here
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

        MUi1k = gener_discrete_dist(semisize, p);
        CPT = zeros(1,2*semisize);
        CPT(1:semisize) = 1-MUi1k;
        CPT(semisize+1:2*semisize) = MUi1k;
        bnet_miss.CPD{N+order(i)} = tabular_CPD (bnet_miss, N+order(i), CPT);

  end

end
