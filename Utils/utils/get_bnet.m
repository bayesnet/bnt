function bnet = get_bnet(data,dag,cpd_name)

% this function can get a bnet from data and dag
%
% bnet = get_bnet(data,dag,cpd_name)
%
% input:
%     data: a matrix data，just support num
%     dag: a dag
%     cpd_name: CPD，tabular means discrete distribution and gaussian means
%               gaussian distribution.
% output:
%     bnet
%
% if you find bug, you can send email to me.
%
% Written by WANGXin(growlithe1205@gmail.com)
%
% 这个函数可以从数据和有向无环图中获得一个网络。
%
% 输入:
%    data: 一个矩阵数据，只支持数字；
%    dag: 一个有向无环图；
%    cpd_name: 条件概率分布方式，tabular代表离散分布 或 gaussian代表高斯分布。
% 输出:
%    bnet
%
% 如果你发现了bug，你可以给我发email
%
% Written by WANGXin(growlithe1205@gmail.com)

[row,~] = size(data);
node_sizes = get_node_sizes(data);
bnet = mk_bnet(dag, node_sizes);

switch cpd_name
    case 'tabular'
        for i=1:row
            bnet.CPD{i} = tabular_CPD(bnet, i);
        end
    case 'gaussian'
        for i=1:row
            bnet.CPD{i} = gaussian_CPD(bnet, i);
        end
end

end