% Test the code using the dag in Fig 1 of Jensen, Jensen, Dittmer, 
% "From influence diagrams to junction trees", UAI 94

% By reverse enginering Fig 2, we infer that the following arcs should
% be absent from the original dag:  b->d1, e->d2, f->d2, g->d4
a=1; b=2; d1=3; c=4; d=5; e=6; f=7; g=8; d2=9; d4=10; i=11; h=12; d3=13; l=14; j=15; k=16;
dag=zeros(16);
dag(a,c)=1;
%dag(b,[c d d1])=1;
dag(b,[c d])=1;
dag(d1,d)=1;
dag(c,e)=1;
dag(d,[e f])=1;
%dag(e,[g d2])=1;
dag(e,[g])=1;
%dag(f,[d2 h])=1;
dag(f,[h])=1;
%dag(g,[d4 i])=1;
dag(g,[i])=1;
dag(d2,i)=1;
dag(d4,l)=1;
dag(i,l)=1;
dag(h,[j k])=1;
dag(d3,k)=1;


[MG, moral_edges] = moralize(dag);
MG(j,k)=1; MG(k,j)=1;  % simulate having a common utility child
% MG now equals fig 2
order = [l j k i h a c d d4 g d3 d2 f e d1 b];
[MTG, cliques, fill_ins] = triangulate(MG, order);
% MTG equals fig 3
ns = 2*ones(1,16);
[jtree, root, cliques2] = mk_strong_jtree(cliques, ns, order, MTG);
jtree2 = mk_rooted_tree(jtree, root);
% jtree2 equals fig 4, with their arrows reversed
