function [bnet,cpdag,BIC_score,nloop] = learn_struct_ges_EM(bnet, data, max_loop, loop_em)
% GES_EM algorithm, learn Bayesian network equivalence classes from incomplete data.
% [bnet,cpdag,BIC_score,nloop] = learn_struct_ges_EM(bnet, data, max_loop)
%
% 2/11/2006: hanene.borchani@gmail.com, francois.olivier.c.h@gmail.com

[N ncases] = size(data);
ns = bnet.node_sizes;
cpdag = dag_to_cpdag(bnet.dag);
order=zeros(1,N);
order=1:N;
tiny = exp(-700);
improve_factor = 0.001; %stop search if the current expected BIC score is less than old_score+old_score*improve_factor,  
log_value = log(ncases); 
update_samples = cell(N, ncases); 
loop = 0;
converged = 0;
if nargin<4, loop_em=6; end

% First phase : INSERT
while ~converged & loop < max_loop
    loop=loop+1;

    fprintf('\n Loop INSERT number %d \n\n', loop);
    engine = jtree_inf_engine(bnet);
    [bnet, LOGLIKE] = learn_params_em(engine, data, loop_em);

    for i=1:N
      s = struct(bnet.CPD{i});
      counts = s.counts(:);
      ll(i) = sum(log(s.CPT(:) + tiny) .* counts);
    end  
    [D,d] = compute_bnet_nparams(bnet);

    [pdags,nodes] =mk_nbrs_of_pdag_add(cpdag,engine); %generate the neighbors of the current network 
    nbrs = length(pdags);
    if nbrs==0
       break;
    end
    dags = pdag_to_dag(pdags); %extract the set of consistent extensions
    nbrs_dags = length(dags);

    if nbrs_dags==0
        disp('The number of consistent extensions of the pdags neighbors is equal to 0');
        break;
    else
      c=0;
      for k=1:nbrs
        if ~isempty(dags{k})
           c=c+1;
           cdags{c}=dags{k};
           cnodes(c,:) = nodes(k,:);
        end
      end
      if c==0
        fprintf('There is no consistent extension');  
        break;
      end
    end

    [ess,LL] = compute_ess(bnet,data,cdags); %compute the estimations for all consistent extensions based on the current bnet
    bic_score0 = sum(LL);
    bic_score0 = bic_score0 - 0.5 * D * log_value; %expected BIC score of the current bnet 
    fprintf('The expected BIC score of the current bnet, bic_score0 =%8.3f \n', bic_score0); 
    bic_score = -inf*ones(1,c); 

    for i=1:c  %compute the expected BIC score of each consistent extension neighbor
        for compt=1:N
            new_CPT{i,compt}=[];
        end
        new_LL = LL;
        new_d=d;
             for num=1:N
                  if ~isempty(ess{i,num}) %consider only nodes whose parent set has been changed 
                       par= mysetdiff(ess{i,num}.domain, num);
                       new_d(num)= prod([ns(par) ns(num)-1]); %compute the new dimension of each node
                       approx = ess{i,num}.counts; 
                       indnum = find(ess{i,num}.domain > num);
                       ldom=length(ess{i,num}.domain);
                       lindnum = length(indnum);
                       if lindnum>0
                          approx = permute(approx, [1:ldom-lindnum-1, indnum(1):ldom, indnum(1)-1]);
                       end
                       new_CPT{i,num} = mk_stochastic(approx);
                       new_LL(num) = sum(log(new_CPT{i,num}(:) + tiny) .* approx(:)); %compute the new LL of each node
                   end
              end
         new_D = sum(new_d);
         bic_score(i) = sum(new_LL) - 0.5 * new_D * log_value; %deduce the expected BIC score of each neighbor
         [a, j] = max(bic_score);
     end

     [BIC_score, best] = max(bic_score);
     fprintf('End computing of the expected BIC scores of all neighbors, the maximal one is Bic_score =%8.3f  \n', BIC_score);
     temp = abs(bic_score0) * improve_factor;   %search will finish when the improvment of the expected BIC score is less than 0.1% compare with the previous best result

   if BIC_score > (bic_score0 + temp)
      best_dag = cdags{best};
      new_order=topological_sort(best_dag);
      for j = 1:N
      row = new_order(j);
         for k = 1:ncases
           update_samples{j,k} = data{row,k};
         end
      end
      data = update_samples;

      forbidden=[];
      reversed =[];
      for j=1:N 
            old_parents= sort((find(bnet.dag(:,j)==1))');
            new_parents= sort((find(best_dag(:,j)==1))');
            if ~isequal(new_parents,old_parents)
               reversed =[reversed,j];
            end
      end

      new_dag = best_dag(new_order, new_order); %reshape the best DAG according to new_order
      ns = ns(new_order); %modify the order of ns
      CPDs = bnet.CPD;
      bnet = mk_bnet(new_dag, ns); %make the new best BN structure

      %randomly set the CPD values of each node
      for j=1:N   
        bnet.CPD{j} = tabular_CPD(bnet, j, 'prior_type', 'dirichlet', 'dirichlet_weight', 0);
      end

      lreversed= length(reversed);
      if lreversed ~=0
         for r=1:lreversed  %update the CPDs of nodes whose parent set has been altered using the saved new_CPT
            reverse=find(new_order==reversed(r));
            forbidden=[forbidden,reverse];
            bnet.CPD{reverse} = set_fields(bnet.CPD{reverse}, 'CPT', new_CPT{best,reversed(r)});
         end
      end

      bnet.CPD = copy_CPD(bnet.CPD, CPDs, new_order, forbidden); %copy the CPDs of remaining nodes from CPDs

      cpdag = dag_to_cpdag(bnet.dag); %get the best equivalence class

      clear bic_score D d new_D new_d;    %new computations for each iteration
      clear ess LL new_LL new_CPT;
      clear pdags dags cdags nodes cnodes;
   else
      fprintf('No improvement of the expected bic score: End of add phase \n \n'); 
      BIC_score = bic_score0;   
      converged=1;
   end
end

nloop=loop;

loop=0; 
converged=0;

if sum(sum(cpdag))==0, converged=1; else cpdag, end

% Second phase : Delete
 fprintf('Start of delete phase \n ');
 while ~converged & loop < max_loop
    loop=loop+1;
    fprintf('\n Loop DELETE number %d \n\n', loop);
    engine = jtree_inf_engine(bnet);
    [bnet, LOGLIKE] = learn_params_em(engine, data, loop_em);  
    for i=1:N
      s = struct(bnet.CPD{i});
      counts = s.counts(:);
      ll(i) = sum(log(s.CPT(:) + tiny) .* counts);
    end  

    [D,d] = compute_bnet_nparams(bnet);

    [pdags,nodes] = mk_nbrs_of_pdag_del(cpdag,engine);  %generate the neighbors of the current network
    nbrs = length(pdags);
    if nbrs==0
        disp('The number of pdag neighbors is equal to 0');
        break;
    end

    dags = pdag_to_dag(pdags); %extract the set of consistent extensions
    nbrs_dags = length(dags);
    if nbrs_dags==0
        disp('The number of consistent extensions of the pdag neighbors is equal to 0');
        break;
    else
      c=0;
      for k=1:nbrs
        if ~isempty(dags{k})
           c=c+1;
           cdags{c}=dags{k};
           cnodes(c,:) = nodes(k,:);
        end
      end
      if c==0
        fprintf('There is no consistent extension');  
        return;
      end
    end

    [ess,LL] = compute_ess(bnet, data,cdags); %compute the estimations for all consistent extensions based on the current bnet
    bic_score0 = sum(LL);
    bic_score0 = bic_score0 - 0.5 * D * log_value;  %expected BIC score of the current bnet  
    fprintf('The expected BIC score of the current bnet, bic_score0 =%8.3f \n', bic_score0); 
    bic_score = -inf*ones(1,c); 

    for i=1:c  %Compute the Expected BIC score for each consistent extension neighbor
        for compt=1:N
           new_CPT{i,compt}=[];
        end
        new_LL = LL;
        new_d=d;
             for num=1:N
                  if ~isempty(ess{i,num}) %consider only nodes whose parent set has been changed 
                       par= mysetdiff(ess{i,num}.domain, num);
                       new_d(num)= prod([ns(par) ns(num)-1]); %compute the new dimension of each node
                       approx = ess{i,num}.counts;
                       indnum = find(ess{i,num}.domain > num);
                       ldom=length(ess{i,num}.domain);
                       lindnum = length(indnum);
                       if lindnum>0
                          approx = permute(approx, [1:ldom-lindnum-1, indnum(1):ldom, indnum(1)-1]);
                       end
                       new_CPT{i,num} = mk_stochastic(approx);
                       new_LL(num) = sum(log(new_CPT{i,num}(:) + tiny) .* approx(:)); %compute the new LL of each node
                   end
              end
         new_D = sum(new_d);
         bic_score(i) = sum(new_LL) - 0.5 * new_D * log_value; %deduce the expected BIC score of each neighbor
     end

     [BIC_score, best] = max(bic_score);
     fprintf('End computing of the expected BIC scores of all neighbor, the maximal one is Bic_score =%8.3f  \n', BIC_score);
     temp = abs(bic_score0) * improve_factor;   

   if BIC_score > (bic_score0 + temp)
      best_dag = cdags{best};
      new_order=topological_sort(best_dag);  
      for j = 1:N
      row = new_order(j);
         for k = 1:ncases
           update_samples{j,k} = data{row,k};
         end
      end
      data = update_samples;

      forbidden=[];
      reversed =[];
      for j=1:N 
            old_parents= sort((find(bnet.dag(:,j)==1))');
            new_parents= sort((find(best_dag(:,j)==1))');
            if ~isequal(new_parents,old_parents)
               reversed =[reversed,j];
            end
      end

      new_dag = best_dag(new_order, new_order); %reshape the best DAG according to new_order
      ns = ns(new_order); %modify the order of ns
      CPDs = bnet.CPD;

      bnet = mk_bnet(new_dag, ns); %make the new best BN structure
      %randomly set the CPD values of each node
      for j=1:N   
        bnet.CPD{j} = tabular_CPD(bnet, j, 'prior_type', 'dirichlet', 'dirichlet_weight', 0);
      end

      lreversed= length(reversed);
      if lreversed ~=0
        for r=1:lreversed  %update the CPDs of nodes whose parent set has been altered using the saved new_CPT
            reverse=find(new_order==reversed(r));
            forbidden=[forbidden,reverse];
            bnet.CPD{reverse} = set_fields(bnet.CPD{reverse}, 'CPT', new_CPT{best,reversed(r)});
         end
      end

      bnet.CPD = copy_CPD(bnet.CPD, CPDs, new_order, forbidden); %copy the CPDs of remaining nodes from CPDs

      cpdag = dag_to_cpdag(bnet.dag); %get the best equivalence class

      clear bic_score D d new_D new_d;    % new computations for each iteration
      clear ess LL new_LL new_CPT;
      clear pdags dags cdags nodes cnodes;
   else
      fprintf('No improvement of the expected bic score : End of delete phase \n \n'); 
      BIC_score = bic_score0;   
      converged=1;
   end
end

cpdag

nloop=nloop+loop; %total iteration number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D,d]= compute_bnet_nparams(bnet)
N = length(bnet.dag);
d = zeros(1,N);
  for i=1:N
    a = struct(bnet.CPD{i});
    d(i) = a.nparams;
  end
D = sum(d);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ess, LL] = compute_ess(bnet, samplesM, dags)
%compute the estimations for all consistent extensions based on the current bnet
tiny = exp(-700);
N = length(bnet.dag);
ns = bnet.node_sizes;
ncases = size(samplesM, 2);
ess0 = cell(1,N);
LL = zeros(1, N); 
nbrs= length(dags); 
copy=zeros(nbrs,N);

for i =1:nbrs
   for x=1:N
       ess{i,x}=[];
   end
end

for i =1:N
   parents = bnet.parents{i};
   family = [parents, i];
   ess0{i} = 0 * myones(ns(family));
end

for i =1:nbrs
     neighbor=dags{i};
     for j=1:N 
         old_parents = sort((find(bnet.dag(:,j)==1))');
         new_parents = sort((find(neighbor(:,j)==1))');
         if ~isequal(new_parents,old_parents)
             domain = sort(myunion(new_parents, j));
             copy(i,j)= find_same_domain(ess,j,domain,i);
             ess{i,j}.domain = domain;
             ess{i,j}.counts = 0 * myones(ns(domain));
         end
     end
end 

engine = jtree_inf_engine(bnet);

for l =1:ncases
    evidence = samplesM(:, l);
    [engine, ll] = enter_evidence(engine, evidence);
    ns_eff = ns;
    ns_eff(~isemptycell(evidence)) = 1;
    Vmarg = cell(1,N);
    for i =1:N
       Vmarg{i} = marginal_nodes(engine, i);
    end
    for i = 1:N
       parents = bnet.parents{i};
       family = [parents, i];
       nfamily = length(family);
       Fmarg = [];
       for j = 1:nfamily
          Fmarg = multiply_one_marginal(Fmarg, Vmarg{family(j)}, ns_eff);
       end
       fullm = add_ev_to_dmarginal(Fmarg, evidence, ns);
       ess0{i} = ess0{i} + fullm.T;
    end

  for i = 1:nbrs
       for j=1:N
          if ~isempty(ess{i,j})
              if ~copy(i,j)
                 domain = ess{i,j}.domain;
                 Fmarg = []; 
                 for compteur=1:length(domain)
                      Fmarg = multiply_one_marginal(Fmarg, Vmarg{domain(compteur)}, ns_eff);
                 end
                 fullm = add_ev_to_dmarginal(Fmarg, evidence, ns);
                 ess{i,j}.counts = ess{i,j}.counts + fullm.T; 
              end
          end
       end
  end 
  clear Vmarg;
end 

for i = 1:nbrs
    for j=1:N
      if copy(i,j)
       ess{i,j}.counts = ess{copy(i,j),j}.counts;
      end
    end
end

for i=1:N
   s = struct(bnet.CPD{i});
   counts = ess0{i};
   LL(i) = sum(log(s.CPT(:) + tiny) .* counts(:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = find_same_domain(ess, j, domain, length)
index = 0;
for i = 1:length
   if ~isempty(ess{i,j})
       if isequal(domain, ess{i,j}.domain)
          index= i;
          break;
       end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  newCPD = copy_CPD(newCPD, CPDs, order, forbidden)
%copy CPDs from old bnet to best bnet, except those nodes whose parent set has been altered
N = length(order);
for i=1:N
    if ~mysubset(i,forbidden)
      a = order(i);
      s = struct(CPDs{a});
      CPT = s.CPT;
      newCPD{i} = set_fields(newCPD{i}, 'CPT', CPT);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
