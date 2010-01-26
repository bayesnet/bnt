function [Gs, op, nodes] = mk_nbrs_of_digraph(G0)
% MK_NBRS_OF_DIGRAPH Make all digraphs that differ from G0 by a single edge deletion, addition or reversal
% [Gs, op, nodes] = mk_nbrs_of_digraph(G0)
%
% Gs(:,:,i) is the i'th neighbor
% op{i} = 'add', 'del', or 'rev' is the operation used to create the i'th neighbor. 
% nodes(i,1:2) are the head and tail of the operated-on arc.

debug = 0; % the vectorized version is about 3 to 10 times faster

n = length(G0);
[I,J] = find(G0); % I(k), J(k) is the k'th edge
E = length(I); % num edges present in G0

% SINGLE EDGE DELETIONS

Grep = repmat(G0(:), 1, E); % each column is a copy of G0
% edge_ndx(k) is the scalar location of the k'th edge 
edge_ndx = find(G0);
% edge_ndx = subv2ind([n n], [I J]); % equivalent
% We set (ndx(k), k) to 0 for k=1:E in Grep
ndx = subv2ind(size(Grep), [edge_ndx(:) (1:E)']);
G1 = Grep;
G1(ndx) = 0;
Gdel = reshape(G1, [n n E]);


% if debug
% % Non-vectorized version
% ctr = 1;
% for e=1:E
%   i = I(e); j = J(e);
%   Gdel2(:,:,ctr) = G0;
%   Gdel2(i,j,ctr) = 0;
%   ctr = ctr + 1;
% end
% assert(isequal(Gdel, Gdel2));
% end


% SINGLE EDGE REVERSALS

% rev_edge_ndx(k) is the scalar location of the k'th reversed edge
%rev_edge_ndx = find(G0'); % different order to edge_ndx, which is bad
rev_edge_ndx = subv2ind([n n], [J I]);
% We set (rev_edge_ndx(k), k) to 1 for k=1:E in G1
% We have already deleted i->j in the previous step
ndx = subv2ind(size(Grep), [rev_edge_ndx(:) (1:E)']);
G1(ndx) = 1;
Grev = reshape(G1, [n n E]);

% if debug
% % Non-vectorized version
% ctr = 1;
% for e=1:E
%   i = I(e); j = J(e);
%   Grev2(:,:,ctr) = G0;
%   Grev2(i,j,ctr) = 0;
%   Grev2(j,i,ctr) = 1;
%   ctr = ctr + 1;
% end
% assert(isequal(Grev, Grev2));
% end


% SINGLE EDGE ADDITIONS

Gbar = ~G0; % Gbar(i,j)=1 iff there is no i->j edge in G0
Gbar = setdiag(Gbar, 0); % turn off self loops
[Ibar,Jbar] = find(Gbar); 

bar_edge_ndx = find(Gbar);
Ebar = length(Ibar); % num edges present in Gbar
Grep = repmat(G0(:), 1, Ebar); % each column is a copy of G0
ndx = subv2ind(size(Grep), [bar_edge_ndx(:) (1:Ebar)']);
Grep(ndx) = 1;
Gadd = reshape(Grep, [n n Ebar]);

% if debug
% % Non-vectorized version
% ctr = 1;
% for e=1:length(Ibar)
%   i = Ibar(e); j = Jbar(e);
%   Gadd2(:,:,ctr) = G0;
%   Gadd2(i,j,ctr) = 1;
%   ctr = ctr + 1;
% end
% assert(isequal(Gadd, Gadd2));
% end


Gs = cat(3, Gdel, Grev, Gadd);

nodes = [I J;
	 I J;
	 Ibar Jbar];

op = cell(1, E+E+Ebar);
op(1:E) = {'del'};
op(E+1:2*E) = {'rev'};
op(2*E+1:end) = {'add'};


% numeric output:
% op(i) = 1, 2, or 3, if the i'th neighbor was created by adding, deleting or reversing an arc.

ADD = 1;
DEL = 2;
REV = 3;

%op = [repmat(DEL, 1, E) repmat(REV, 1, E) repmat(ADD, 1, Ebar)];
