function [bnet, order, BIC_score, LOGLIKE] = learn_struct_EM(bnet, samplesM, max_loop)
% LEARN_STRUCT_EM(), structural EM algorithm , learn structure and parameters
% from missing data.
% [bnet, order, BIC_score] = learn_struct_EM(bnet, samplesM, max_loop)

tiny = exp(-700);
improve_factor = 0.001;     %when current BIC score is less than old_score+old_score*improve_factor, stop search 
N = length(bnet.dag); 
ncases = size(samplesM, 2);
log_value = log(ncases);
ns = bnet.node_sizes;
dag = zeros(N,N);
order = zeros(1,N);    %save the label of each node in current dag correspond to the original dag
order = 1:N;           %original dag has nodes label 1:N
CPT = cell(2);         %save the modified node's CPT of the bnet that has the highest score in an iteration
                       %in cases "del" and "add", there is only one CPT will be save, in "rev" need to save two CPTs.
update_samples = cell(N, ncases);   %in this algorithm, because the label of the next dag will be different with the
                                    %last dag, so the label of the trainning data will be modified, too

loop = 0;
evidence = cell(1,N);
while loop<max_loop    % generally set the max_loop to 30   
   loop = loop + 1
   engine = jtree_inf_engine(bnet);
   [bnet, LOGLIKE] = learn_params_em(engine, samplesM, 10);     % default set the parameter EM runs 10 iterations
   for i=1:N
      s = struct(bnet.CPD{i});
      counts = s.counts(:);
      ll(i) = sum(log(s.CPT(:) + tiny) .* counts);
   end   
   [D,d] = compute_bnet_nparams(bnet);
   
   [nbrs, ops, nodes, orders] = mk_nbrs_of_dag_topo(bnet.dag); 
   nGs = length(nbrs);
   
   [ec, ec1, LL] = compute_approx_ess(bnet, samplesM, ops, nodes);
   bic_score0 = sum(LL);
   bic_score0 = bic_score0 - 0.5 * D * log_value;  % bic score of current bnet
 
   bic_score = zeros(1,nGs);   % save each neighbour dag(bnet)'s bic score 
   for i=1:nGs
      bic_score(i) = -inf;
   end
   for i=1:nGs
      edge = nodes(i,:);
      switch ops{i}
      case 'del'
         head = edge(1);
         tail = edge(2);
         approx_ess = ec{i}.counts;
         CPT1 = mk_stochastic(approx_ess); 
         
         LL1 = LL;
         LL1(tail) = sum(log(CPT1(:) + tiny) .* approx_ess(:));
         d1 = d;
         d1(tail) = d(tail) / ns(head);
         D1 = sum(d1);
         bic_score(i) = sum(LL1) - 0.5 * D1 * log_value;
         [a, j] = max(bic_score);
         if j==i                                     % if the current dag has the highest bic score, save it's CPT(s)
            CPT{1} = CPT1;
         end

      case 'add'
         head = edge(1);
         tail = edge(2);
         approx_ess = ec{i}.counts;
         if head>tail                     % now, the "ess" is in ascent manner, accord with the labels in the "domain" field.
            n = length(ec{i}.domain);     % need permute , so that "ess" contain the last dimension is about the "tail" node.
            approx_ess = permute(approx_ess, [1:n-2, n, n-1]);       % because there is only one "edge" modified, only need 
         end                                                         % to exchange the last two dimension if needed.
         CPT1 = mk_stochastic(approx_ess);

         LL1 = LL;
         d1 = d;
         LL1(tail) = sum(log(CPT1(:) + tiny) .* approx_ess(:));
         d1(tail) = d(tail) * ns(head);
         D1 = sum(d1);
         bic_score(i) = sum(LL1) - 0.5 * D1 * log_value;
         [a, j] = max(bic_score);
         if j==i
            CPT{1} = CPT1;
         end

      case 'rev'          % ops "rev" influent two family, equals the combination of a "del" and an "add"
         % "del" an edge
         head = edge(1);
         tail = edge(2);
         approx_ess = ec1{i}.counts;
         CPT1 = mk_stochastic(approx_ess); 
         LL1 = LL;
         LL1(tail) = sum(log(CPT1(:) + tiny) .* approx_ess(:));
         d1 = d;
         d1(tail) = d(tail) / ns(head);

         % "add" an edge
         head = edge(2);
         tail = edge(1);
         approx_ess = ec{i}.counts;
         if head>tail                     % now, the "ess" is in ascent manner, accord with the labels in the "domain" field.
            n = length(ec{i}.domain);     % need permute , so that "ess" contain the last dimension is about the "tail" node.
            approx_ess = permute(approx_ess, [1:n-2, n, n-1]);       % because there is only one "edge" modified, only need 
         end                                                         % to exchange the last two dimension if needed.
         CPT2 = mk_stochastic(approx_ess);
         LL1(tail) = sum(log(CPT2(:) + tiny) .* approx_ess(:));
         d1(tail) = d(tail) * ns(head);

         D1 = sum(d1);
         bic_score(i) = sum(LL1) - 0.5 * D1 * log_value;
         [a, j] = max(bic_score);
         if j==i
            CPT{1} = CPT1;
            CPT{2} = CPT2;
         end
      end
   end

   [BIC_score, i] = max(bic_score);
   temp = abs(bic_score0) * improve_factor;      % search will be finish when the improvment of bic score 
                                        % less than 0.1% compare with the previous best result
   if BIC_score > (bic_score0 + temp)
      dag1 = nbrs{i};                   % new best dag
      order1 = orders{i};               % labels of each nodes altered from last iteration

      % the labels of each nodes are altered, so the "data" will need to "re-arrange" according to the new order
      for j = 1:N
         row = order1(j);
         for k = 1:ncases
            update_samples{j,k} = samplesM{row,k};
         end
      end
      samplesM = update_samples;

      dag = dag1(order1, order1);       % "reshape" the best dag, make it as an "upper trianglar"
      ns = ns(order1);                  % also must modify the order of "ns"
      CPDs = bnet.CPD;
      bnet = mk_bnet(dag, ns);          % use the best dag now to produce a new bnet, with altered nodes labels
      for j=1:N                         % randomly set the CPTs values of each CPDs
         bnet.CPD{j} = tabular_CPD(bnet, j, 'prior_type', 'dirichlet', 'dirichlet_weight', 0);
      end
      edge = nodes(i,:);

      % update the CPDs of new best bnet(dag) using corresponding CPDs of last iteration.
      % copy the old CPTs that not altered. 
      % set the altered CPTs from the saved "CPT" variables 
      switch ops{i}
      case 'del'
         tail = edge(2);
         tail = find(order1==tail);
         bnet.CPD{tail} = set_fields(bnet.CPD{tail}, 'CPT', CPT{1});
         forbidden = [tail];
         bnet.CPD = copy_CPD(bnet.CPD, CPDs, order1, forbidden);
      case 'add'
         tail = edge(2);
         tail = find(order1==tail);
         bnet.CPD{tail} = set_fields(bnet.CPD{tail}, 'CPT', CPT{1});
         forbidden = [tail];
         bnet.CPD = copy_CPD(bnet.CPD, CPDs, order1, forbidden);
      case 'rev'
         head = edge(2);
         head = find(order1==head);
         bnet.CPD{head} = set_fields(bnet.CPD{head}, 'CPT', CPT{1});
         tail = edge(1);
         tail = find(order1==tail);
         bnet.CPD{tail} = set_fields(bnet.CPD{tail}, 'CPT', CPT{2});
         forbidden = [head, tail];
         bnet.CPD = copy_CPD(bnet.CPD, CPDs, order1, forbidden);
      end

      % draw a graph for the new best dag with nodes labels are the same as the original
      order = order(order1);
%      labels = cellstr(int2str(order'));
%      figure(loop+1);
%      draw_graph(dag,labels);

      clear bic_score D d;           % for each iteration, re-compute all the expected counts and bic score
      clear ec ec1 LL;
      clear nbrs ops nodes orders;
   else
      BIC_score = bic_score0;    % if there is no improvement in bic score, stop the search, and return
      break;
   end
end
BIC_score


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D,d]= compute_bnet_nparams(bnet)
%
%
N = length(bnet.dag);
d = zeros(1,N);
for i=1:N
   a = struct(bnet.CPD{i});
   d(i) = a.nparams;
end
D = sum(d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  newCPD = copy_CPD(newCPD, CPDs, order, forbidden)
%copy CPDs from old bnet to new bnet, except those nodes has been modified
%
N = length(order);
for i=1:N
   if ~mysubset(i, forbidden)
      a = order(i);
      s = struct(CPDs{a});
      CPT = s.CPT;
      newCPD{i} = set_fields(newCPD{i}, 'CPT', CPT);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ec, ec1, LL] = compute_approx_ess(bnet, samplesM, ops, nodes)
%compute all neighbours' needed approximate ess based on current bnet.
%
tiny = exp(-700);
N = length(bnet.dag);
ns = bnet.node_sizes;
ncases = size(samplesM, 2);
nGs = length(ops);
ec0 = cell(1,N);
ec = cell(1, nGs);         %store each neighbours' altered family's approximate ess.
ec1 = cell(1, nGs);        %since operator 'rev' need to alter two families, ec1 store the ess of family deleted an edge
copy = zeros(1, nGs);
copy1 = zeros(1, nGs);
LL = zeros(1, N);          %For current bnet, LL store each nodes's LogLike based on approximate ess.
for i =1:nGs
   ec{i}.domain = [];
   ec{i}.counts = [];
   ec1{i}.domain = [];
   ec1{i}.counts = [];
end
for i =1:N
   parents = bnet.parents{i};
   family = [parents, i];
   ec0{i} = 0 * myones(ns(family));
end
for i =1:nGs
   edge = nodes(i, :);
   switch ops{i}
   case 'del'
      head = edge(1);
      tail = edge(2);
      parents = bnet.parents{tail};
      parents = mysetdiff(parents, head);
      domain = [parents, tail];
      copy(i) = find_same_domain(ec, domain, i);
      ec{i}.domain = domain;
      ec{i}.counts = 0 * myones(ns(domain));
   case 'add'
      head = edge(1);
      tail = edge(2);
      parents = bnet.parents{tail};
      parents = [parents, head, tail];
      domain = sort(parents);
      copy(i) = find_same_domain(ec, domain, i);
      ec{i}.domain = domain;
      ec{i}.counts = 0 * myones(ns(domain));
   case 'rev'
      head = edge(1);
      tail = edge(2);
      parents = bnet.parents{tail};
      parents = mysetdiff(parents, head);
      domain = [parents, tail];
      copy1(i) = find_same_domain(ec, domain, i);
      ec1{i}.domain = domain;
      ec1{i}.counts = 0 * myones(ns(domain));

      head = edge(2);
      tail = edge(1);
      parents = bnet.parents{tail};
      parents = [parents, head, tail];
      domain = sort(parents);
      copy(i) = find_same_domain(ec, domain, i);
      ec{i}.domain = domain;
      ec{i}.counts = 0 * myones(ns(domain));
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
      ec0{i} = ec0{i} + fullm.T;
   end

   for i = 1:nGs
      switch ops{i}
      case 'del'
         if ~copy(i)
            domain = ec{i}.domain;
            Fmarg = [];
            for j=1:length(domain)
               Fmarg = multiply_one_marginal(Fmarg, Vmarg{domain(j)}, ns_eff);
            end
            fullm = add_ev_to_dmarginal(Fmarg, evidence, ns);
            ec{i}.counts = ec{i}.counts + fullm.T;
         end
      case 'add'
         if ~copy(i) 
            domain = ec{i}.domain;
            Fmarg = [];
            for j=1:length(domain)
               Fmarg = multiply_one_marginal(Fmarg, Vmarg{domain(j)}, ns_eff);
            end
            fullm = add_ev_to_dmarginal(Fmarg, evidence, ns);
            ec{i}.counts = ec{i}.counts + fullm.T;
         end
      case 'rev'
         if ~copy1(i) 
            domain = ec1{i}.domain;
            Fmarg = [];
            for j=1:length(domain)
               Fmarg = multiply_one_marginal(Fmarg, Vmarg{domain(j)}, ns_eff);
            end
            fullm = add_ev_to_dmarginal(Fmarg, evidence, ns);
            ec1{i}.counts = ec1{i}.counts + fullm.T;
         end

         if ~copy(i) 
            domain = ec{i}.domain;
            Fmarg = [];
            for j=1:length(domain)
               Fmarg = multiply_one_marginal(Fmarg, Vmarg{domain(j)}, ns_eff);
            end
            fullm = add_ev_to_dmarginal(Fmarg, evidence, ns);
            ec{i}.counts = ec{i}.counts + fullm.T;
         end
      end
   end
   clear Vmarg;
end

for i =1:nGs
   if copy(i)
      ec{i}.counts = ec{copy(i)}.counts;
   end
   if copy1(i)
      ec1{i}.counts = ec{copy1(i)}.counts;
   end
end

for i=1:N
   s = struct(bnet.CPD{i});
   counts = ec0{i};
   LL(i) = sum(log(s.CPT(:) + tiny) .* counts(:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = find_same_domain(ec, domain, length)
%
%
index = 0;
for i = 1:length
   if isequal(domain, ec{i}.domain)
      index = i;
      break;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





