function  [cpdag] = dag_to_cpdag1(dags)
% 2
% (also works with a cell array of dags, returning a cell array of cpdags)
% DAG_TO_CPDAG produce a N*N matrix which values respect :
%
% 	If the edge is compelled then 1 on the edge.
%	If the edge is reversible then 1 on the edge and 1 in the reverse edge.
%
% Make sure that the entry is a DAG.
%
% See D.M. Chickering: "Learning Equivalence Classes of Bayesian-Network Structures".
%
% 
% francois.olivier.c.h@gmail.com, philippe.leray@univ-nantes.fr

if ~iscell(dags)
    dag=cell(1,1);
    dag{1}=dags;
else
    dag=dags;
end

for da=1:length(dag)
    cpdags{da} = abs(label_edges(dag{da}));
end

if ~iscell(dags)
    cpdag=cpdags{1};
else
    cpdag=cpdags;
end

%%==============================================================================

function [label] = label_edges(dag)
% LABEL-EDGES produce a N*N matrix which values are
% 	+1 if the edge is compelled or
%	-1 if the edge is reversible.
% Make sure that the entry is a DAG.
%
% francois.olivier.c.h@gmail.com

N=length(dag);
[order xedge yedge] = order_edges(dag);
label = 2*dag; % all edges as unknown

NbEdges = length(xedge) ;
%xedge=x, yedge=y, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for Edge=1:NbEdges,
    xlow=xedge(Edge);
    ylow=yedge(Edge);
    if label(xlow,ylow)==2
        fin = 0;

        %wcompelled = find(label(:,xlow)~=2);
        wcompelled = find(label(:,xlow)==1);

        parenty = find(label(:,ylow)~=0);
        %sonsy = find(label(ylow,:)~=0);

        for s = 1:length(wcompelled)
           w = wcompelled(s);
           if ~ismember(w,parenty)
                label(xlow,ylow)=1; 
                label(ylow,xlow)=0; %
                label(parenty,ylow)=1;
                label(ylow,parenty)=0; %
                %label(ylow,sonsy)=1; %
                %label(sonsy,ylow)=0; %
                fin = 1;

            elseif fin == 0
                label(w,ylow)=1;
                %label(ylow,w)=0; %
            end
        end
        if fin == 0
            parentx = [xlow ; find(label(:,xlow)~=0)];
            if ~isempty(mysetdiff(parenty,parentx))
                %label(xlow,ylow)=1; %
                %label(ylow,xlow)=0; %

                label(find(label(:,ylow)==2),ylow)=1;
                label(ylow,find(label(ylow,:)==2))=1; %
            else	
                label(xlow,ylow)=-1;
                label(ylow,xlow)=-1; %
                ttp = find(label(:,ylow)==2);
                label(ttp,ylow)=-1;
                label(ylow,ttp)=-1; %
            end
        end
    end
end

%%%========================================================================================
function [order, x, y] = order_edges(dag)
% ORDER_EDGES produce a total (natural) ordering over the edges in a DAG.
% Make sure that the entry is a DAG.
%
% francois.olivier.c.h@gmail.com
%
% 2 mai 2003

if acyclic(dag)==0
    error('Requires an acyclic graph');
end

N=length(dag);
order = zeros(N,N);

node_order = topological_sort(dag);
[tmp oo] = sort(node_order);

dag2=dag(node_order,node_order);
[x y]=find(flipud(dag2)==1);
nb_edges=length(x);

if nb_edges~=0
  order(sub2ind([N N],N+1-x,y))=1:nb_edges ;
end

order=order(oo,oo);
x=node_order(N+1-x);
y=node_order(y);

