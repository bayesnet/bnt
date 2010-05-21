function resu = ismemberclique(v,cliques)
% b = ismemberclique(v,cliques)
%

finiclique = 0 ; resu=0 ; 
cl=1;  ncl=length(cliques) ;

while (~finiclique) & (cl<=ncl);
    if ismember(v,cliques{cl})
        resu=1;
        finiclique=1 ;
    end
    cl=cl+1;
end
