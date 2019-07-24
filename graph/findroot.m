function root = findroot(bnet, cliques)

%% findroot is to find the strong root in a clique tree assume it has one
%% in the tree. For a clique tree constructed from a strongly triangulated
%% graph, an interface clique that contains all discrete parents
%% and at least one continuous node from a connected continuous component
%% is for sure to be available as a guaranteed strong root.
%% -By Wei Sun, George Mason University, 4/17/2010.

%% We choose the interface clique that contains the max number 
%% of interface nodes to be the strong root.
n0 = 0 ;
for i=1:length(cliques)
    % check hybrid cliques
    hc = intersect(cliques{i}, bnet.cnodes) ; 
    hd = intersect(cliques{i}, bnet.dnodes) ;
    if ~isempty(hd) & ~isempty(hc)
        nd = length(hd) ;
        if nd > n0
            root = i ;
            n0 = nd ;
        end
    end    
end
