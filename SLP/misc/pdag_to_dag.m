function G2 = pdag_to_dag(pdags)
% (also works with a cell array of pdags, returning a cell array of dags)
% dag = pdag_to_dag(pdag)
%
% cf Dor and Tarsi (1992) :
%    A simple algorithm to construct a consistent extention of a partially oriented graph.
%
% francois.olivier.c.h@gmail.com


if ~iscell(pdags)
    pdag=cell(1,1);
    pdag{1}=pdags;
else
    pdag=pdags;
end

for da=1:length(pdag)
    %fprintf('%d ',da)
    G=pdag{da};
    G2 = G; A = G;
    N = size(G,1);
    empty_loop = 0;

    while ~isempty(find(A))
        [x x_y_undirected] = select_vertex(A);

        if x==0
            fprintf('pdag_to_dag error : This pdag does not admit any extension.\n');
            G2=[];
            break
        end
        G2(x,x_y_undirected) = 0; G2(x_y_undirected,x) = 1;

        A(x,:) = 0;
        A(:,x) = 0;
    end
    dags{da}=G2;
end
if ~iscell(pdags)
    G2=dags{1};
else
    G2=dags;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol, x_y] = select_vertex(G)
N = size(G,1);
sol = 0;
x = 0;
fini=0 ;
while ~fini

    x = x+1;
    if x>N
        fini=1;
    else
        beforex=find(G(:,x));
        afterx=find(G(x,:));

        if ~(isempty(beforex)&isempty(afterx))
            x_y = myintersect(beforex,afterx);
            beforex=mysetdiff(beforex,x_y);
            afterx=mysetdiff(afterx,x_y);

            if isempty(afterx) % x is a sink
                Ax=  myunion(x_y,beforex);
                for y=x_y
                    % Adjacents of y
                    Ay = myunion(find(G(:,y)), find(G(y,:)));
                    Ay = myunion(Ay,y);
                    if isempty(setdiff(Ax,Ay))
                        fini=fini+1; 
                    else
                        break;
                    end
                end
                if fini==length(x_y)
                    sol=x; fini=1;
                else
                    fini=0;
                end
            end
        end
    end
end % while
