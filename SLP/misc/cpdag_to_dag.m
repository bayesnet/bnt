function dag = cpdag_to_dag(cpdags)
% dags = cpdag_to_dag(cpdags)
%
% CPDAG_TO_DAG produce a N*N matrix of a dag which instantiate cpdag.
% (also works with a cell array of cpdags, returning a cell array of dags)
% make sur that your entry is a completed PDAG
% this function can't be use instead of PDAG_TO_DAG
%
%  see Chickering (2002) : Learning equivalence classes of bayesian networks, JMLR2, pp475-479
%
% francois.olivier.c.h@gmail.com, philippe.leray@univ-nantes.fr
% 31 march 2006

if ~iscell(cpdags)
    cpdag=cell(1,1);
    cpdag{1}=cpdags;
else
    cpdag=cpdags;
end

for da=1:length(cpdag)

    N=length(cpdag{da});
    dag=cpdag{da};

    unprocessed = find_nodes_in_undirected_component(dag);
    while ~isempty(unprocessed)
      nbr_parents = [];
      for i=1:length(unprocessed)
        nbr_parents(end+1)=length(parents(dag-dag.*dag',unprocessed(i))); %nbr_parents(end+1)=length(parents(dag,unprocessed(i)));
      end
      [tmp, idx] = max(nbr_parents);
      node = unprocessed(idx);
      dag(parents(dag.*dag',node),node)=0; %dag(parents(dag,node),node)=0;
      %dag(myintersect(parents(dag,node),unprocessed), node)=0; % Wei Lu
      unprocessed = find_nodes_in_undirected_component(dag);
    end

    dags{da}=dag;
end

if ~iscell(cpdags)
    dag=dags{1};
else
    dag=dags;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function unprocessed = find_nodes_in_undirected_component(dag)
undirected_edges = dag.*dag';
[unprocessed, tmp] = find(undirected_edges);
unprocessed = unique(unprocessed);
