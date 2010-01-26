% example to illustrate why nodes must be numbered topologically.
% Due to Shinya OHTANI <ohtani@pdp.crl.sony.co.jp>
% 9 June 2004

%%%%%%%%% WRONG RESULTS because 2 -> 1 
% should have P(parent|no evidence) = prior = [03. 0.7]

node = struct('ChildNode', 1, ...
               'ParentNode', 2);

adjacency = zeros(2);
adjacency([node.ParentNode], node.ChildNode) = 1;

value = {{'TRUE'; 'FALSE'}, ...
          {'TRUE'; 'FALSE'}};

bnet = mk_bnet(adjacency, [2 2]);
bnet.CPD{node.ChildNode} = tabular_CPD(bnet, node.ChildNode, [0.2 0.4 0.8 0.6]);
bnet.CPD{node.ParentNode} = tabular_CPD(bnet, node.ParentNode, [0.3 0.7]);

evidence = cell(1,2);
% evidence{node.ChildNode} = 1;
% evidence{node.ParentNode} = 1;

engine = jtree_inf_engine(bnet);
[engine, loglik] = enter_evidence(engine, evidence);


marg = marginal_nodes(engine, node.ChildNode);
disp(sprintf('    ChildNode     : %8.6f   %8.6f',marg.T(1),marg.T(2)) );
marg = marginal_nodes(engine, node.ParentNode);
disp(sprintf('    ParentNode    : %8.6f   %8.6f',marg.T(1),marg.T(2)) );

% 
%     ChildNode   : 0.534483   0.465517
%     ParentNode  : 0.155172   0.844828
% loglik = 0.15



%%%%%%%%% RIGHT RESULTS because 1 -> 2

node = struct('ChildNode', 2, ...
               'ParentNode', 1);


adjacency = zeros(2);
adjacency([node.ParentNode], node.ChildNode) = 1;

value = {{'TRUE'; 'FALSE'}, ...
          {'TRUE'; 'FALSE'}};

bnet = mk_bnet(adjacency, [2 2]);
bnet.CPD{node.ChildNode} = tabular_CPD(bnet, node.ChildNode, [0.2 0.4 0.8 0.6]);
bnet.CPD{node.ParentNode} = tabular_CPD(bnet, node.ParentNode, [0.3 0.7]);

evidence = cell(1,2);
% evidence{node.ChildNode} = 1;
% evidence{node.ParentNode} = 1;

engine = jtree_inf_engine(bnet);
[engine, loglik] = enter_evidence(engine, evidence);


marg = marginal_nodes(engine, node.ChildNode);
disp(sprintf('    ChildNode     : %8.6f   %8.6f',marg.T(1),marg.T(2)) );
marg = marginal_nodes(engine, node.ParentNode);
disp(sprintf('    ParentNode    : %8.6f   %8.6f',marg.T(1),marg.T(2)) );
