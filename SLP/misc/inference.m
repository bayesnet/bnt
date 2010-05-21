function [proba_post, engin2]= inference(bnet, data, node)
% Make bayesian inference on data
% [proba_post, engine]= inference(bnet, data, node)
%
% INPUTS :
% - bnet, the structure of the bayesian network gived by mk_bnet.
% - data(i,m), node i in case m.
% - node, the node we interrogating.
%
% OUTPUTS :
% - proba_post, the posteriors probabilities.
% - engine, the inference engine.
%
% francois.olivier.c.h@gmail.com

engine=jtree_inf_engine(bnet);
[N L]=size(data);
proba_post=zeros(L,bnet.node_sizes(node));
for i=1:L
   evidence(1:N)=data(1:N,i);
   evidence{node}=[];
   [engin2, ll]=enter_evidence(engine,evidence);
   marg=marginal_nodes(engin2,node);
   proba_post(i,:)=marg.T';
end