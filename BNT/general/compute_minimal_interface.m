function clqs = compute_minimal_interface(intra, inter)

int = compute_fwd_interface(intra, inter);
ss  = length(intra);
Z = zeros(ss);
dag = [intra inter;
       Z     intra];
G = moralize(dag);
intra2 = G(1:ss,1:ss);
inter2 = G(1:ss,(1:ss)+ss);
G = unroll_dbn_topology(intra2, inter2, ss);
T = ss;
last_slice = (1:ss) + (T-1)*ss;
G = (G + G')/2; % mk symmetric
G2 = (expm(full(G)) > 0); % closure of graph
G3 = G2(last_slice, last_slice);
[c,v] = scc(G3); % connected components
ncomp = size(v,1);
clqs = cell(1,ncomp);
for i=1:ncomp
  ndx = find(v(i,:)>0);
  clqs{i} = v(i,ndx);
  clqs{i} = myintersect(clqs{i}, int);
end

