
function [eg] = dag_to_essential_graph(dag)
cpdag = dag_to_cpdag(dag);
eg = dag + dag .* (cpdag + cpdag');

return;




% Coverts a DAG into Essential Graph where edges are coded by 2 and 3, 2 is
% directed edge and 3 is bidirected edge and is at one (the same as the original DAG) of the two
% symetrical places. 

% Is implemented by the algorithm of Max Chickering in D.M.Chickering (1995). 
% A transformational characterization of equivalent Bayesian network structures. 
% In Proceedings of Eleventh Conference on Uncertainty in Artificial Intelligence, Montreal, QU,
% pages 87-98. Morgan Kaufmann 
% http://research.microsoft.com/~dmax/publications/uai95.pdf 

% Implemented by Tomas Kocka, AAU.

function [eg] = dag_to_essential_graph(dagx)

%print_dag(dagx); % Just checking input

order = topological_sort(dagx); % get the topological order of nodes and their number

% fprintf('the topological order is: %d',order);
% fprintf('\n');

[nx,ny] = size(dagx); % gets the number of nodes, note that nx == ny
[I,J] = find(dagx); % finds all nonzero elements in the adjacency matrix, i.e. arcs in the DAG - however we will overwrite it in a special order
% we will sort the arcs from lowest possible y and highest possible x, arcs are x->y
e = 1;
for y = 1:ny
    for x = nx:-1:1
        %fprintf('x %d ',order(x)); fprintf('y %d ',order(y));
        if dagx(order(x),order(y)) == 1 
            I(e) = order(x);
            J(e) = order(y);
            e = e + 1;
            %fprintf('x order %d',x);
            %fprintf('y order %d',y);
            %fprintf('\n');
        end
    end
end


% fprintf('the arcs are: %d',I);
% fprintf('\n');
% fprintf('the arcs are: %d',J);
% fprintf('\n');


% Now we have to decide which arcs are part of the essential graph and
% which are undirected edges in the essential graph.
% Undecided arc in the DAG are 1, directed in EG are 2 and undirected in EG
% are 3.


for e = 1:length(I)
    if dagx(I(e),J(e)) == 1
        cont = true;
        for w = 1:nx 
            if dagx(w,I(e)) == 2
                if dagx(w,J(e)) ~= 0
                    dagx(w,J(e)) = 2;
                else
                    for ww = 1:nx
                        if dagx(ww,J(e)) ~= 0
                           dagx(ww,J(e)) = 2;
                        end
                    end % and now skip the rest and start with another arc from the list
                    w = nx;
                    cont = false;
                end
            end
        end
        if cont
           exists = false;
           for z = 1:nx
               %fprintf('test %d',dagx(z,J(e)));
               if dagx(z,J(e)) ~= 0 & z ~= I(e) & dagx(z,I(e)) == 0
                  exists = true; 
                  for ww = 1:nx
                        if dagx(ww,J(e)) == 1
                           dagx(ww,J(e)) = 2;
                        end 
                  end
               end
           end
           if ~ exists
               for ww = 1:nx
                   if dagx(ww,J(e)) == 1
                      dagx(ww,J(e)) = 3;
                   end 
               end  
           end
        end
    end            
end

%print_dag(dagx); % Just checking output







