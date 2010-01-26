function bnet = mk_markov_chain_bnet(N, Q)

dag = zeros(N);
dag(1,2)=1; dag(2,3)=1;
ns = Q*ones(1,N); 
bnet = mk_bnet(dag, ns);
for i=1:N
  bnet.CPD{i} = tabular_CPD(bnet, i);
end
