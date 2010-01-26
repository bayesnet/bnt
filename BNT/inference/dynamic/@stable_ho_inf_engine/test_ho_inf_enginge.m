function [engine,engine2] = test_ho_inf_enginge(order,T)

assert(order >= 1)
% Model a SISO system, i. e. all node are one-dimensional
% The nodes are numbered as follows
% u(t) = 1 input
% y(t) = 2 model output
% z(t) = 3 noise
% q(t) = 4 observed output = noise + model output

ns = [1 1 1 1];

% Model a linear system, i.e. there are no discrete nodes
dn = [];

% Modeling of connections within a time slice
intra = zeros(4);
intra(2,4) = 1; % Connection y(t) -> q(t)
intra(3,4) = 1; % Connection z(t) -> q(t)

% Connections to the next time slice
inter = zeros(4,4,order);
inter(1,2,1) = 1; % u(t) -> y(t+1);
inter(2,2,1) = 1; %y(t) -> y(t+1);
inter(3,3,1) = 1; %z(t) -> z(t+1);

if order >= 2
    inter(1,2,2) = 1; % u(t) -> y(t+2);
    inter(2,2,2) = 1; % y(t) -> y(t+2);
end

for i = 3: order
    inter(:,:,i) = inter(:,:,i-1); %u(t) -> y(t+i) y(t) -> y(t) +i
end;


% Compution of a higer order Markov Model
bnet = mk_higher_order_dbn(intra,inter,ns,'discrete',dn);
bnet2 = mk_dbn(intra,inter(:,:,1),ns,'discrete',dn)


%Calculation of the number of nodes with different parameters
%There is one input and one output nodes  2
%There are two different disturbance node 2
%There are order +1 nodes for y           1 + order
numOfNodes = 5 + order; 

% First input node
bnet.CPD{1} = gaussian_CPD(bnet,1,'mean',0);
bnet2.CPD{1} = gaussian_CPD(bnet,1,'mean',0);
% Modeled output
bnet.CPD{2} = gaussian_CPD(bnet,2,'mean',0);
bnet2.CPD{2} = gaussian_CPD(bnet,2,'mean',0);
%Disturbance
bnet.CPD{3} = gaussian_CPD(bnet,3,'mean',0);
bnet2.CPD{3} = gaussian_CPD(bnet,3,'mean',0);

%Qutput
bnet.CPD{4} = gaussian_CPD(bnet,4,'mean',0);
bnet2.CPD{4} = gaussian_CPD(bnet,4,'mean',0);


%Output node in the second time-slice
%Remember that node number 6 is an example for 
%the fifth equivalence class
bnet.CPD{5} = gaussian_CPD(bnet,6,'mean',0);
bnet2.CPD{5} = gaussian_CPD(bnet,6,'mean',0);

%Disturbance node in the second time slice
bnet.CPD{6} = gaussian_CPD(bnet,7,'mean',0);
bnet2.CPD{6} = gaussian_CPD(bnet,7,'mean',0);

% Modeling of the remaining nodes for y
for i = 7:numOfNodes
    bnet.CPD{i} = gaussian_CPD(bnet,(i - 6)*4 + 7,'mean',0);
end

% Generation of the inference engine
engine = dv_unrolled_dbn_inf_engine(bnet,T);
engine2 = jtree_unrolled_dbn_inf_engine(bnet,T);







