function [bnet, onodes]=hme_topobuilder(nodes_info);
%
% HME topology builder
%
% ----------------------------------------------------------------------------------------------------
% -> pierpaolo_b@hotmail.com   or   -> pampo@interfree.it
% ----------------------------------------------------------------------------------------------------

nodes_num=size(nodes_info,2);
dag = zeros(nodes_num);
list=[1:nodes_num];
for i=1:(nodes_num-1)
    app=[];
    app=list((i+1):end);
    dag(i,app) = 1;
end
onodes = [1 nodes_num];                           
dnodes = list(2:end-1);
if nodes_info(1,end)>0,
    dnodes=[dnodes nodes_num];
end
ns = nodes_info(2,:);

bnet = mk_bnet(dag, ns, dnodes);
clamped = 0;

bnet.CPD{1} = root_CPD(bnet, 1);

rand('state', 50);
randn('state', 50);

for i=2:nodes_num,
    if (nodes_info(1,i)==0)&(nodes_info(4,i)==1),
        bnet.CPD{i} = gaussian_CPD(bnet, i, [], [], [], 'full');
    elseif (nodes_info(1,i)==0)&(nodes_info(4,i)==2),
        bnet.CPD{i} = gaussian_CPD(bnet, i, [], [], [], 'diag');
    elseif (nodes_info(1,i)==0)&(nodes_info(4,i)==3),
        bnet.CPD{i} = gaussian_CPD(bnet, i, [], [], [], 'full', 'tied');
    elseif (nodes_info(1,i)==0)&(nodes_info(4,i)==4),
        bnet.CPD{i} = gaussian_CPD(bnet, i, [], [], [], 'diag', 'tied');        
    elseif nodes_info(1,i)==1,
        %bnet.CPD{i} = dsoftmax_CPD(bnet, i, [], [], clamped, nodes_info(4,i));
	bnet.CPD{i} = softmax_CPD(bnet, i, 'clamped', clamped, 'max_iter', nodes_info(4,i));
    else
        bnet.CPD{i} = mlp_CPD(bnet, i, nodes_info(3,i), [], [], [], [], clamped, nodes_info(4,i));
    end
end
