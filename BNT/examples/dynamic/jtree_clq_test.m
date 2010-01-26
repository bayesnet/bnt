% Construct various DBNs and examine their clique structure.
% This was used to generate various figures in chap 3-4 of my thesis.

% Examine the cliques in the unrolled mildew net

%dbn = mk_mildew_dbn;
dbn = mk_chmm(4);
ss = dbn.nnodes_per_slice;
T = 7;
N = ss*T;
bnet = dbn_to_bnet(dbn, T);

constrained = 0;
if constrained
  stages = num2cell(unroll_set(1:ss, ss, T), 1);
else
  stages = { 1:N; };
end
clusters = {};
%[jtree, root, cliques, B, w, elim_order, moral_edges, fill_in_edges] = ...
%    dag_to_jtree(bnet, bnet.observed, stages, clusters);
[jtree, root, cliques] =  graph_to_jtree(moralize(bnet.dag), ones(1,N), stages, clusters);

flip=1;
clf;[dummyx, dummyy, h] = draw_dbn(dbn.intra, dbn.inter, flip, T, -1);
dir = '/home/eecs/murphyk/WP/Thesis/Figures/Inf/MildewUnrolled';
mk_ps_from_clqs(dbn, T, cliques, [])
%mk_collage_from_clqs(dir, cliques)


% Examine the cliques in the cascade DBN

% A-A
%  \
% B B
%  \
% C C
%  \
% D D
ss = 4;
intra = zeros(ss);
inter = zeros(ss);
inter(1, [1 2])=1;
for i=2:ss-1
  inter(i,i+1)=1;
end


% 2 coupled HMMs 1,3  and 2,4
ss = 4;
intra = zeros(ss);
inter = zeros(ss); % no persistent edges
%inter = diag(ones(ss,1)); % persitence edges
inter(1,3)=1; inter(3,1)=1;
inter(2,4)=1; inter(4,2)=1;

%bnet = mk_fhmm(3);
bnet = mk_chmm(4);
intra = bnet.intra;
inter = bnet.inter;

clqs = compute_minimal_interface(intra, inter);
celldisp(clqs)




% A A
%  \
% B B
%  \
% C C
%  \
% D-D
ss = 4;
intra = zeros(ss);
inter = zeros(ss);
for i=1:ss-1
  inter(i,i+1)=1;
end
inter(4,4)=1;



ns = 2*ones(1,ss);
dbn = mk_dbn(intra, inter, ns);
for i=2*ss
  dbn.CPD{i} = tabular_CPD(bnet, i);
end

T = 4;
N = ss*T;
bnet = dbn_to_bnet(dbn, T);

constrained = 1;
if constrained
  % elim first 3 slices first in any order
  stages = {1:12, 13:16};
  %stages = num2cell(unroll_set(1:ss, ss, T), 1);
else
  stages = { 1:N; };
end
clusters = {};
%[jtree, root, cliques, B, w, elim_order, moral_edges, fill_in_edges] = ...
%    dag_to_jtree(bnet, bnet.observed, stages, clusters);
[jtree, root, cliques] =  graph_to_jtree(moralize(bnet.dag), ones(1,N), stages, clusters);





% Examine the cliques in the 1.5 slice DBN

%dbn = mk_mildew_dbn;
dbn = mk_water_dbn;
%dbn = mk_bat_dbn;
ss = dbn.nnodes_per_slice;
int = compute_fwd_interface(dbn);
bnet15 = mk_slice_and_half_dbn(dbn, int);
N = length(bnet15.dag);
stages = {1:N};

% bat
%cl1 = [16 17 19 7 14];
%cl2 = [27 25 21 23 20];
%clusters = {cl1, cl2, cl1+ss, cl2+ss};

% water
%cl1 = 1:2; cl2 = 3:6; cl3 = 7:8;
%clusters = {cl1, cl2, cl3, cl1+ss, cl2+ss, cl3+ss};

%clusters = {};
clusters = {int, int+ss};
%[jtree, root, cliques, B, w, elim_order, moral_edges, fill_in_edges] = ...
%    dag_to_jtree(bnet15, bnet.observed, stages, clusters);
[jtree, root, cliques] =  graph_to_jtree(moralize(bnet15.dag), ones(1,N), stages, clusters);

clq_len = [];
for c=1:length(cliques)
  clq_len(c) = length(cliques{c});
end
hist(clq_len, 1:max(clq_len));
h=hist(clq_len, 1:max(clq_len));
axis([1 max(clq_len)+1 0 max(h)+1])
xlabel('clique size','fontsize',16)
ylabel('number','fontsize',16)




