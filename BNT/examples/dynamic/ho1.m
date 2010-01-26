function ho1()

% Example of how to create a higher order DBN
% Written by Rainer Deventer <deventer@informatik.uni-erlangen.de> 3/28/03

bnet = createBNetNL();

%%%%%%%%%%%%


function bnet = createBNetNL(varargin)
     % Generate a Bayesian network, which is able to model nonlinearities at
% the input. The only input is the order of the dynamic system. If this 
% parameter is missing, the an order of two is assumed
if nargin > 0 
    order = varargin{1}
else
    order = 2;
end

ss = 6; % For each time slice the following nodes are modeled
        % ud(t_k) Discrete node, which decides whether saturation is reached.
        %         Node number 2
        % uv(t_k) Visible input node with node number  2
        % uh(t_k) Hidden  input node with node number 3     
        % y(t_k)  Modeled output, Number 4
        % z(t_k)  Disturbing variable, number 5
        % q(t_k), number6 6

intra = zeros(ss,ss);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Within each timeslice ud(t_k) is connected with uv(t_k) and uh(t_k)    %
% This part is used to model saturation                                  %
% A connection from  uv(t_k) to uh(t_k) is omitted                       %
% Additionally   y(t_k) is connected with q(t_k). To model the disturbing%
% value z(t_k) is connected with q(t_k).                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intra(1,2:3) = 1; % Connections ud(t_k) -> uv(t_k) and ud(t_k) -> uh(t_k)
intra(4:5,6) = 1; % Connectios  y(t_k)  -> q(t_k)  and z(t_k)  -> q(t_k) 


  
inter = zeros(ss,ss,order);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Markov assumption is not met as connections from time slice t to t+2 %
% exist.                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:order
    if i == 1
        inter(1,1,i) = 1; %Connect the discrete nodes. This is necessary to improve
                          %the disturbing reaction
        inter(3,4,i) = 1; %Connect uh(t_{k-1}) with y(t_k)
        inter(4,4,i) = 1; %Connect y(t_{k-1})  with y(t_k)    
        inter(5,5,i) = 1; %Connect z(t_{k-1})  with z(t_k)
    else
        inter(3,4,i) = 1; %Connect uh(t_{k-i}) with y(t_k)
        inter(4,4,i) = 1; %Connect  y(t_{k-i}) with y(t_k)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the dimensions of the discrete nodes. Node 1 has two states     %
% 1 = lower saturation reached                                           %
% 2 = Upper saturation reached                                           %
% Values in between are model by probabilities between 0 and 1           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
node_sizes = ones(1,ss);
node_sizes(1) = 2;
dnodes = [1];

eclass = [1:6;7 2:3 8 9 6;7 2:3 10 11 6];
bnet = mk_higher_order_dbn(intra,inter,node_sizes,...
                           'discrete',dnodes,...
                           'eclass',eclass);

cov_high = 400;
cov_low  = 0.01;
weight1 = randn(1,1);
weight2 = randn(1,1);
weight3 = randn(1,1);
weight4 = randn(1,1);

numOfNodes = 5 + order;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nodes of the first time-slice   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete input node, 
bnet.CPD{1} = tabular_CPD(bnet,1,'CPT',[1/2 1/2],'adjustable',0);


% Modeled visible input
bnet.CPD{2} = gaussian_CPD(bnet,2,'mean',[0 10],'clamp_mean',1,...
                            'cov',[10 10],'clamp_cov',1);

% Modeled hidden input
bnet.CPD{3} = gaussian_CPD(bnet,3,'mean',[0, 10],'clamp_mean',1,...
			          'cov',[0.1 0.1],'clamp_cov',1);

% Modeled output in the first timeslice, thus there are no parents
% Usuallz the output nodes get a low covariance. But in the first
% time-slice a prediction of the output is not possible due to 
% missing information
bnet.CPD{4} = gaussian_CPD(bnet,4,'mean',0,'clamp_mean',1,...
			          'cov',cov_high,'clamp_cov',1);

%Disturbance
bnet.CPD{5} = gaussian_CPD(bnet,5,'mean',0,...
                                  'cov',[4],...
                                  'clamp_mean',1,...
                                  'clamp_cov',1);

%Observed output. 
bnet.CPD{6} = gaussian_CPD(bnet,6,'mean',0,...
                                  'clamp_mean',1,...
                                  'cov',cov_low,'clamp_cov',1,...
                                  'weights',[1 1],'clamp_weights',1);

% Discrete node at second time slice
bnet.CPD{7} = tabular_CPD(bnet,7,'CPT',[0.6 0.4 0.4 0.6],'adjustable',0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Node for the model output %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnet.CPD{8} = gaussian_CPD(bnet,10,'mean',0,...
				   'cov',cov_high,...
				   'clamp_mean',1,...
			           'clamp_cov',1);
%                                   'weights',[0.0791 0.9578]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Node for the disturbance %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
bnet.CPD{9} = gaussian_CPD(bnet,11,'mean',0,'clamp_mean',1,...
                                   'cov',[4],'clamp_cov',1,...
                                   'weights',[1],'clamp_weights',1);
                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Node for the model output %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnet.CPD{10} = gaussian_CPD(bnet,16,'mean',0,'clamp_mean',1,...
                                    'cov',cov_low,'clamp_cov',1);
%                                   'weights',[0.0188 -0.0067 0.0791 0.9578]);



    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Node for the disturbance %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
bnet.CPD{11} = gaussian_CPD(bnet,17,'mean',0,'clamp_mean',1,...
                                           'cov',[0.2],'clamp_cov',1,...
                                           'weights',[1],'clamp_weights',1);




