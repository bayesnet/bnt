function [new_nbrs, new_ops, new_nodes, new_topos] = mk_nbrs_of_dag_topo(G0)
% MK_NBRS_OF_DAG_TOPO Make all DAGs that differ from G0 by a single edge deletion, addition or reversal
% [new_nbrs, new_ops, new_nodes, new_topos] = mk_nbrs_of_dag_topo(G0)
%
% new_nbrs{i} is the i'th neighbor of G0.
% new_ops{i} = 'add', 'del', or 'rev' is the operation used to create the i'th neighbor.
% new_nodes(i,1:2) are the head and tail of the operated-on arc.
% new_topos are topological orders of the neighbours.
%
% We implement the fast acyclicity check described by P. Giudici and R. Castelo,
% "Improving MCMC model search for data mining", submitted to J. Machine Learning, 2001.
%
% Written by Qian Diao <qian.diao@intel.com> on 19 Nov 01
% Reference are ..\BNT\graph\mk_nbrs_of_dag.m, ..\BNT\learning\learn_struct_mcmc.m and 
% ..\BNT\graph\topological_sort.m
% Copyright Intel 2001
%
new_nbrs = {};
new_ops = {};
new_nodes = [];
new_topos = {};
cs = {};

n = length(G0);
indeg = zeros(1,n);
zero_indeg = []; % a stack of nodes with no parents
for i=1:n
  indeg(i) = length(parents(G0,i));
  cs{i} = children(G0, i); 
  if indeg(i)==0
    zero_indeg = [i zero_indeg];
  end
end

dag = G0;
[nbrs, ops, nodes] = mk_nbrs_of_digraph(dag);
A = init_ancestor_matrix(dag);
%assert(acyclic(new_dag));

d1 = 1;
for d = 1:length(ops)
  i = nodes(d, 1); j = nodes(d, 2);
  legal = 0;
  switch ops{d}
   case 'add',
    if A(i,j)==0
      legal = 1;
    end
   case 'del',
    legal = 1;
   case 'rev',
    ps = mysetdiff(parents(dag, j), i);
    % if any(A(ps,i)) then there is a path i -> parent of j -> j
    % so reversing i->j would create a cycle
    legal = ~any(A(ps, i));
  end

  if legal 
    tmp = nbrs(:,:,d);
    new_nbrs{d1} = tmp;
    new_ops{d1} =  ops{d};
    new_nodes(d1,1:2) = nodes(d,1:2);

    % obtain the topological orders of neighbour dags
    zero_indeg_nbr = [];
    indeg_nbr = [];
    cs_nbr = [];

    switch ops{d}
     case 'add' % i is a new parent of j 
      zero_indeg_nbr = zero_indeg;
      indeg_nbr = indeg;  
      if ~isempty(find(zero_indeg == j)) % j is not a root anymore   
        zero_indeg_nbr = mysetdiff(zero_indeg, j); 
      end   
      indeg_nbr(j) = indeg(j)+1;

      t_nbr=1;
      order_nbr = zeros(1,n);
      while ~isempty(zero_indeg_nbr)
        v_nbr = zero_indeg_nbr(1); % pop v
        zero_indeg_nbr = zero_indeg_nbr(2:end);
        order_nbr(t_nbr) = v_nbr;
        t_nbr = t_nbr + 1;
        if v_nbr == i % j is a new child of i
          cs_nbr = sort([j cs{i}]);
        else
          cs_nbr = cs{v_nbr};
        end  
        for k = 1:length(cs_nbr)
          c_nbr = cs_nbr(k);
          indeg_nbr(c_nbr) = indeg_nbr(c_nbr) - 1;
          if indeg_nbr(c_nbr) == 0
            zero_indeg_nbr = [c_nbr zero_indeg_nbr]; % push c 
          end
        end
      end

     case 'del' % i is not a parent of j anymore
      zero_indeg_nbr = zero_indeg; 
      indeg_nbr = indeg;  
      if length(parents(tmp, j))==0  
        zero_indeg_nbr = -sort(-[zero_indeg, j]); % descending order
      end   
      indeg_nbr(j) = indeg(j) - 1;

      t_nbr=1;
      order_nbr = zeros(1,n);
      while ~isempty(zero_indeg_nbr)
        v_nbr = zero_indeg_nbr(1); % pop v
        zero_indeg_nbr = zero_indeg_nbr(2:end);
        order_nbr(t_nbr) = v_nbr;
        t_nbr = t_nbr + 1;
        if v_nbr == i % j is not a child of i anymore
          cs_nbr = mysetdiff(cs{i}, j);
        else
          cs_nbr = cs{v_nbr};
        end  
        for k = 1:length(cs_nbr)
          c_nbr = cs_nbr(k);
          indeg_nbr(c_nbr) = indeg_nbr(c_nbr) - 1;
          if indeg_nbr(c_nbr) == 0
            zero_indeg_nbr = [c_nbr zero_indeg_nbr]; % push c 
          end
        end
      end

     case 'rev' %i is a new child of j and j is a new parent of i
      zero_indeg_nbr = zero_indeg;  
      indeg_nbr = indeg;
      if ~isempty(find(zero_indeg == i))    
        zero_indeg_nbr = mysetdiff(zero_indeg_nbr, i); 
      end   
      if length(parents(tmp, j))==0  
        zero_indeg_nbr = -sort(-[zero_indeg_nbr, j]); % decending order
      end   
      indeg_nbr(i) = indeg(i)+1;
      indeg_nbr(j) = indeg(j)-1;

      t_nbr=1;
      order_nbr = zeros(1,n);
      while ~isempty(zero_indeg_nbr)
        v_nbr = zero_indeg_nbr(1); % pop v
        zero_indeg_nbr = zero_indeg_nbr(2:end);
        order_nbr(t_nbr) = v_nbr;
        t_nbr = t_nbr + 1;
        cs_nbr = cs{v_nbr};
        if v_nbr == i % j is not a child of i anymore
          cs_nbr = mysetdiff(cs{i}, j);
        end
        if v_nbr == j % i is a new child of j  
          cs_nbr = sort([i cs{j}]); 
        end  
        for k = 1:length(cs_nbr)
          c_nbr = cs_nbr(k);
          indeg_nbr(c_nbr) = indeg_nbr(c_nbr) - 1;
          if indeg_nbr(c_nbr) == 0
            zero_indeg_nbr = [c_nbr zero_indeg_nbr]; % push c 
          end
        end
      end
     end  

    new_topos{d1} = order_nbr; 
    d1 = d1+1;
  end 
end

clear nbrs ops nodes;



%%%%%%%%%
function A = update_row(A, j, dag)
% We compute row j of A
A(j, :) = 0;
ps = parents(dag, j);
if ~isempty(ps)
  A(j, ps) = 1;
end
for k=ps(:)'
  anck = find(A(k,:));
  if ~isempty(anck)
    A(j, anck) = 1;
  end
end

%%%%%%%%
function A = init_ancestor_matrix(dag)
order = topological_sort(dag);
A = zeros(length(dag));
for j=order(:)'
  A = update_row(A, j, dag);
end

   