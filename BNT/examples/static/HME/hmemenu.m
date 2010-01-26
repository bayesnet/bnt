% dataset      -> (1=>user data) or (2=>toy example)
% type         -> (1=> Regression model) or (2=>Classification model)
% num_glevel   -> number of hidden nodes in the net (gating levels)
% num_exp      -> number of experts in the net
% branch_fact  -> dimension of the hidden nodes in the net
% cov_dim      -> root node dimension
% res_dim      -> output node dimension
% nodes_info   -> 4 x num_glevel+2 matrix that contain all the info about the nodes:
%                 nodes_info(1,:) = nodes type: (0=>gaussian)or(1=>softmax)or(2=>mlp)
%                 nodes_info(2,:) = nodes size: [cov_dim   num_glevel x branch_fact   res_dim]
%                 nodes_info(3,:) = hidden units number (for mlp nodes)
%                                  |- optimizer iteration number (for softmax & mlp CPD)
%                 nodes_info(4,:) =|- covariance type (for gaussian CPD)-> 
%                                  | (1=>Full)or(2=>Diagonal)or(3=>Full&Tied)or(4=>Diagonal&Tied)
% fh1 -> Figure: data & decizion boundaries; fh2 -> confusion matrix; fh3 -> LL trace                                                                       
% test_data    -> test data matrix
% train_data   -> training data matrix
% ntrain       -> size(train_data,2)
% ntest        -> size(test_data,2)
% cases        -> (cell array) training data formatted for the learning engine
% bnet         -> bayesian net before learning
% bnet2        -> bayesian net after learning
% ll           -> log-likelihood before learning
% LL2          -> log-likelihood trace
% onodes       -> obs nodes in bnet & bnet2
% max_em_iter  -> maximum number of interations of the EM algorithm
% train_result -> prediction on the training set (as test_result)
% 
% IMPORTANT: CHECK the loading path (lines 64 & 364)
% ----------------------------------------------------------------------------------------------------
% -> pierpaolo_b@hotmail.com   or   -> pampo@interfree.it
% ----------------------------------------------------------------------------------------------------

error('this no longer works with the latest version of BNT')

clear all;
clc;
disp('---------------------------------------------------');
disp('  Hierarchical Mixtures of Experts models builder ');
disp('---------------------------------------------------');
disp(' ')
disp('   Using this script you can build both an HME model')
disp('as in [Wat94] and [Jor94] i.e. with ''softmax'' gating')
disp('nodes and ''gaussian'' ( for regression ) or ''softmax''') 
disp('( for classification ) expert node, and its variants')
disp('called ''gated nets'' where we use ''mlp'' models in')
disp('place of a number of ''softmax'' ones [Mor98], [Wei95].')
disp('  You can decide to train and test the model on your')
disp('datasets  or  to evaluate its  performance on  a toy')
disp('example.')
disp(' ')
disp('Reference')
disp('[Mor98] P. Moerland (1998):')
disp('        Localized mixtures of experts. (http://www.idiap.ch/~perry/)')
disp('[Jor94] M.I. Jordan, R.A. Jacobs (1994):')
disp('        HME and the EM algorithm. (http://www.cs.berkeley.edu/~jordan/)')
disp('[Wat94] S.R. Waterhouse, A.J. Robinson (1994):') 
disp('        Classification using HME. (http://www.oigeeza.com/steve/)')
disp('[Wei95] A.S. Weigend, M. Mangeas (1995):') 
disp('        Nonlinear gated experts for time series.')
disp(' ')

if 0
disp('(See the figure)')
pause(5);
%%%%%WARNING!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im_path=which('HMEforMatlab.jpg');
fig=imread(im_path, 'jpg');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Units','pixels','MenuBar','none','NumberTitle','off', 'Name', 'HME model');
image(fig); 
axis image;
axis off;
clear fig;
set(gca,'Position',[0 0 1 1])
disp('(Press any key to continue)')
pause
end

clc
disp('---------------------------------------------------');
disp('              Specify the Architecture             ');
disp('---------------------------------------------------');
disp(' ');
disp('What kind of model do you need?')
disp(' ')
disp('1) Regression ')
disp('2) Classification')
disp(' ')
type=input('1 or 2?: ');
if (isempty(type)|(~ismember(type,[1 2]))), error('Invalid value'); end
clc
disp('----------------------------------------------------');
disp('               Specify the Architecture             ');
disp('----------------------------------------------------');
disp(' ')
disp('Now you have to set the number of experts and gating')
disp('levels in the net.  This script builds only balanced')
disp('hierarchy with the same branching factor (>1)at each')
disp('(gating) level. So remember that: ')
disp(' ')
disp('         num_exp = branch_fact^num_glevel           ')
disp(' ')
disp('with branch_fact >=2.')
disp('You can also set to zeros the number of gating level')
disp('in order to obtain a classical GLM model.           ')
disp(' ')
disp('----------------------------------------------------');
disp(' ')
num_glevel=input('Insert the number of gating levels {0,...,20}: ');
if (isempty(num_glevel)|(~ismember(num_glevel,[0:20]))), error('Invalid value'); end
nodes_info=zeros(4,num_glevel+2);
if num_glevel>0, %------------------------------------------------------------------------------------
    for i=2:num_glevel+1,
        clc
        disp('----------------------------------------------------');
        disp('               Specify the Architecture             ');
        disp('----------------------------------------------------');
        disp(' ')   
        disp(['-> Gating network ', num2str(i-1), ' is a: '])
        disp(' ')
        disp('   1) Softmax model');
        disp('   2) Two layer perceptron model')
        disp(' ')
        nodes_info(1,i)=input('1 or 2?: ');
        if (isempty(nodes_info(1,i))|(~ismember(nodes_info(1,i),[1 2]))), error('Invalid value'); end
        disp(' ')
        if nodes_info(1,i)==2,
           nodes_info(3,i)=input('Insert the number of units in the hidden layer: ');
           if (isempty(nodes_info(3,i))|(floor(nodes_info(3,i))~=nodes_info(3,i))|(nodes_info(3,i)<=0)), 
              error(['Invalid value: ', num2str(nodes_info(3,i)), ' is not a positive integer!']);
           end
           disp(' ')
        end
        nodes_info(4,i)=input('Insert the optimizer iteration number: ');
        if (isempty(nodes_info(4,i))|(floor(nodes_info(4,i))~=nodes_info(4,i))|(nodes_info(4,i)<=0)), 
           error(['Invalid value: ', num2str(nodes_info(4,i)), ' is not a positive integer!']);
        end    
    end
    clc
    disp('---------------------------------------------------------');
    disp('                 Specify the Architecture                ');
    disp('---------------------------------------------------------');
    disp(' ')
    disp('Now you have to set the number  of experts in the network');
    disp('The value will be adjusted in order to obtain a hierarchy');
    disp('as said above.')
    disp(' ');    
    num_exp=input(['Insert the approximative number of experts (>=', num2str(2^num_glevel), '): ']);
    if (isempty(num_exp)|(num_exp<=0)|(num_exp<2^num_glevel)), 
        error('Invalid value');
    end
    app1=0; base=2;
    while app1<num_exp,
        app1=base^num_glevel;
        base=base+1;
    end
    app2=(base-2)^num_glevel;
    branch_fact=base-1;
    if app2>=(2^num_glevel)&(abs(app2-num_exp)<abs(app1-num_exp)),
        branch_fact=base-2;
    end
    clear app1 app2 base;
    disp(' ')
    disp(['The effective number of experts in the net is: ', num2str(branch_fact^num_glevel), '.'])
    disp(' ');
else
    clc
    disp('---------------------------------------------------------');
    disp('        Specify the Architecture (GLM model)             ');
    disp('---------------------------------------------------------');
    disp(' ')
end % END of: if num_glevel>0-------------------------------------------------------------------------

if type==2,
    disp(['-> Expert node is a: '])
    disp(' ')
    disp('   1) Softmax model');
    disp('   2) Two layer perceptron model')
    disp(' ')
    nodes_info(1,end)=input('1 or 2?: ');
    if (isempty(nodes_info(1,end))|(~ismember(nodes_info(1,end),[1 2]))), 
        error('Invalid value'); 
    end
    disp(' ')
    if nodes_info(1,end)==2,
       nodes_info(3,end)=input('Insert the number of units in the hidden layer: ');
       if (isempty(nodes_info(3,end))|(floor(nodes_info(3,end))~=nodes_info(3,end))|(nodes_info(3,end)<=0)), 
           error(['Invalid value: ', num2str(nodes_info(3,end)), ' is not a positive integer!']);
       end
       disp(' ')
    end
    nodes_info(4,end)=input('Insert the optimizer iteration number: ');
    if (isempty(nodes_info(4,end))|(floor(nodes_info(4,end))~=nodes_info(4,end))|(nodes_info(4,end)<=0)), 
        error(['Invalid value: ', num2str(nodes_info(4,end)), ' is not a positive integer!']);
    end
elseif type==1,
    disp('What kind of covariance matrix structure do you want?')
    disp(' ')
    disp('   1) Full');
    disp('   2) Diagonal')
    disp('   3) Full & Tied');
    disp('   4) Diagonal & Tied')

    disp(' ')
    nodes_info(4,end)=input('1, 2, 3 or 4?: ');
    if (isempty(nodes_info(4,end))|(~ismember(nodes_info(4,end),[1 2 3 4]))), 
        error('Invalid value'); 
    end  
end
clc
disp('----------------------------------------------------');
disp('                    Specify the Input               ');
disp('----------------------------------------------------');
disp(' ')
disp('Do you want to...')
disp(' ')
disp('1) ...use your own dataset?')
disp('2) ...apply the model on a toy example?')
disp(' ')
dataset=input('1 or 2?: ');
if (isempty(dataset)|(~ismember(dataset,[1 2]))), error('Invalid value'); end
if dataset==1,
    if type==1,
        clc
        disp('-------------------------------------------------------');
        disp('        Specify the Input - Regression problem         ');
        disp('-------------------------------------------------------');
        disp(' ')
        disp('Be sure that each row of your data matrix is an example');
        disp('with the covariate values that precede the respond ones')
        disp(' ')
        disp('-------------------------------------------------------');
        disp(' ')
        cov_dim=input('Insert the covariate space dimension: ');
        if (isempty(cov_dim)|(floor(cov_dim)~=cov_dim)|(cov_dim<=0)), 
          error(['Invalid value: ', num2str(cov_dim), ' is not a positive integer!']);
        end
        disp(' ')
        res_dim=input('Insert the dimension of the respond variable: ');
        if (isempty(res_dim)|(floor(res_dim)~=res_dim)|(res_dim<=0)), 
            error(['Invalid value: ', num2str(res_dim), ' is not a positive integer!']);
        end 
        disp(' ');
    elseif type==2
        clc
        disp('-------------------------------------------------------');
        disp('      Specify the Input - Classification problem       ');
        disp('-------------------------------------------------------');
        disp(' ')
        disp('Be sure that each row of your data matrix is an example');
        disp('with the covariate values that precede the class labels');
        disp('(integer value >=1).                                   ');
        disp(' ')
        disp('-------------------------------------------------------');
        disp(' ')
        cov_dim=input('Insert the covariate space dimension: ');
        if (isempty(cov_dim)|(floor(cov_dim)~=cov_dim)|(cov_dim<=0)), 
          error(['Invalid value: ', num2str(cov_dim), ' is not a positive integer!']);
        end
        disp(' ')
        res_dim=input('Insert the number of classes: ');
        if (isempty(res_dim)|(floor(res_dim)~=res_dim)|(res_dim<=0)), 
          error(['Invalid value: ', num2str(res_dim), ' is not a positive integer!']);
        end        
        disp(' ')               
    end    
    % ------------------------------------------------------------------------------------------------
    % Loading training data --------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------
    train_path=input('Insert the complete (with extension) path of the training data file:\n >> ','s');    
    if isempty(train_path), error('You must specify a data set for training!'); end
    if ~isempty(findstr('.mat',train_path)),
        ap=load(train_path); app=fieldnames(ap); train_data=eval(['ap.', app{1,1}]);
        clear ap app;
    elseif ~isempty(findstr('.txt',train_path)),
        train_data=load(train_path, '-ascii');
    else
        error('Invalid data format: not a .mat or a .txt file')
    end
    if (size(train_data,2)~=cov_dim+res_dim)&(type==1),
        error(['Invalid data matrix size: ', num2str(size(train_data,2)), ' columns rather than ',...
            num2str(cov_dim+res_dim),'!']);
    elseif (size(train_data,2)~=cov_dim+1)&(type==2),
        error(['Invalid data matrix size: ', num2str(size(train_data,2)), ' columns rather than ',...
            num2str(cov_dim+1),'!']);    
    elseif (~isempty(find(ismember(intersect([train_data(:,end)' 1:res_dim],...
            train_data(:,end)'),[1:res_dim])==0)))&(type==2),
        error('Invalid class label');
    end    
    ntrain=size(train_data,1);
    train_d=train_data(:,1:cov_dim);
    if type==2,
        train_t=zeros(ntrain, res_dim);
        for m=1:res_dim,
            train_t((find(train_data(:,end)==m))',m)=1;
        end
    else
        train_t=train_data(:,cov_dim+1:end);
    end        
    disp(' ')
    % ------------------------------------------------------------------------------------------------
    % Loading test data ------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------
    disp('(If you don''t want to specify a test-set press ''return'' only)');
    test_path=input('Insert the complete (with extension) path of the test data file:\n >> ','s');  
    if ~isempty(test_path),
        if ~isempty(findstr('.mat',test_path)),        
            ap=load(test_path); app=fieldnames(ap); test_data=eval(['ap.', app{1,1}]);
            clear ap app;
        elseif ~isempty(findstr('.txt',test_path)),
            test_data=load(test_path, '-ascii');
        else
            error('Invalid data format: not a .mat or a .txt file')
        end
        if (size(test_data,2)~=cov_dim)&(size(test_data,2)~=cov_dim+res_dim)&(type==1),
            error(['Invalid data matrix size: ', num2str(size(test_data,2)), ' columns rather than ',...
                num2str(cov_dim+res_dim), ' or ', num2str(cov_dim), '!']);
        elseif (size(test_data,2)~=cov_dim)&(size(test_data,2)~=cov_dim+1)&(type==2),
            error(['Invalid data matrix size: ', num2str(size(test_data,2)), ' columns rather than ',...
                num2str(cov_dim+1), ' or ', num2str(cov_dim), '!']);
        elseif (~isempty(find(ismember(intersect([test_data(:,end)' 1:res_dim],...
                test_data(:,end)'),[1:res_dim])==0)))&(type==2)&(size(test_data,2)==cov_dim+1),
            error('Invalid class label');
        end
        ntest=size(test_data,1);        
        test_d=test_data(:,1:cov_dim);
        if (type==2)&(size(test_data,2)>cov_dim),
            test_t=zeros(ntest, res_dim);
            for m=1:res_dim,
                test_t((find(test_data(:,end)==m))',m)=1;
            end
        elseif (type==1)&(size(test_data,2)>cov_dim),
            test_t=test_data(:,cov_dim+1:end);
        end
        disp(' ');
    end
else    
    clc
    disp('----------------------------------------------------');
    disp('                  Specify the Input                 ');
    disp('----------------------------------------------------');
    disp(' ')
    ntrain = input('Insert the number of examples in training (<500): ');
    if (isempty(ntrain)|(floor(ntrain)~=ntrain)|(ntrain<=0)|(ntrain>500)), 
          error(['Invalid value: ', num2str(ntrain), ' is not a positive integer <500!']);
    end        
    disp(' ')
    test_path='toy';
    ntest = input('Insert the number of examples in test (<500): ');
    if (isempty(ntest)|(floor(ntest)~=ntest)|(ntest<=0)|(ntest>500)), 
          error(['Invalid value: ', num2str(ntest), ' is not a positive integer <500!']);
    end        

    if type==2,
        cov_dim=2;
        res_dim=3;
        seed = 42;
        [train_d, ntrain1, ntrain2, train_t]=gen_data(ntrain, seed);
        for m=1:ntrain
            q=[]; q = find(train_t(m,:)==1);
            train_data(m,:)=[train_d(m,:) q];
        end
        [test_d, ntest1, ntest2, test_t]=gen_data(ntest);
        for m=1:ntest
            q=[]; q = find(test_t(m,:)==1);
            test_data(m,:)=[test_d(m,:) q];
        end
    else
        cov_dim=1;
        res_dim=1;
        global HOME
        %%%%%WARNING!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load([HOME '/examples/static/Misc/mixexp_data.txt'], '-ascii');
        %%%%%WARNING!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        train_data = mixexp_data(1:ntrain, :);
        train_d=train_data(:,1:cov_dim); train_t=train_data(:,cov_dim+1:end);
        test_data = mixexp_data(ntrain+1:ntrain+ntest, :);
        test_d=test_data(:,1:cov_dim); 
        if size(test_data,2)>cov_dim,
            test_t=test_data(:,cov_dim+1:end);
        end
    end    
end
% Set the nodes dimension-----------------------------------
if num_glevel>0,
    nodes_info(2,2:num_glevel+1)=branch_fact;
end
nodes_info(2,1)=cov_dim; nodes_info(2,end)=res_dim;
%-----------------------------------------------------------
% Prepare the training data for the learning engine---------
%-----------------------------------------------------------
cases = cell(size(nodes_info,2), ntrain);
for m=1:ntrain,
    cases{1,m}=train_data(m,1:cov_dim)';
    cases{end,m}=train_data(m,cov_dim+1:end)';
end
%-----------------------------------------------------------------------------------------------------
[bnet onodes]=hme_topobuilder(nodes_info);
engine = jtree_inf_engine(bnet, onodes);
clc
disp('---------------------------------------------------------------------');
disp('                         L  E  A  R  N  I  N  G                      ');
disp('---------------------------------------------------------------------');
disp(' ')
ll = 0;
for l=1:ntrain
  scritta=['example number: ', int2str(l),'---------------------------------------------'];
  disp(scritta);
  ev = cases(:,l);
  [engine, loglik] = enter_evidence(engine, ev);
  ll = ll + loglik;
end
disp(' ')
disp(['Log-likelihood before learning: ', num2str(ll)]);
disp(' ')
disp('(Press any key to continue)');
pause
%-----------------------------------------------------------
clc
disp('---------------------------------------------------------------------');
disp('                         L  E  A  R  N  I  N  G                      ');
disp('---------------------------------------------------------------------');
disp(' ')
max_em_iter=input('Insert the maximum number of the EM algorithm iterations: ');
if (isempty(max_em_iter)|(floor(max_em_iter)~=max_em_iter)|(max_em_iter<=1)), 
          error(['Invalid value: ', num2str(ntest), ' is not a positive integer >1!']);
end 
disp(' ')
disp(['Log-likelihood before learning: ', num2str(ll)]);
disp(' ')

[bnet2, LL2] = learn_params_em(engine, cases, max_em_iter);
disp(' ')
fprintf('HME: loglik before learning %f, after %d iters %f\n', ll, length(LL2),  LL2(end));
disp(' ')
disp('(Press any key to continue)');
pause
%-----------------------------------------------------------------------------------
% Classification problem: plot data & decision boundaries if the input data size = 2
% Regression problem: plot data & prediction if the input data size = 1
%-----------------------------------------------------------------------------------
if (type==2)&(nodes_info(2,1)==2)&(~isempty(test_path)),
    fh1=hme_class_plot(bnet2, nodes_info, train_data, test_data);
    disp(' ');
    disp('(See the figure)');
elseif (type==2)&(nodes_info(2,1)==2)&(isempty(test_path)),
    fh1=hme_class_plot(bnet2, nodes_info, train_data);
    disp(' ');
    disp('(See the figure)');
elseif (type==1)&(nodes_info(2,1)==1)&(~isempty(test_path)),
    fh1=hme_reg_plot(bnet2, nodes_info, train_data, test_data);
    disp(' ');
    disp('(See the figure)');
elseif (type==1)&(nodes_info(2,1)==1)&(isempty(test_path)),
    fh1=hme_reg_plot(bnet2, nodes_info, train_data);
    disp(' ')
    disp('(See the figure)');
end
%-----------------------------------------------------------------------------------
% Classification problem: plot confusion matrix
%-----------------------------------------------------------------------------------
if (type==2)
    ztrain=fhme(bnet2, nodes_info, train_d, size(train_d,1));  
    [Htrain, trainRate]=confmat(ztrain, train_t); % CM on the training set
    fh2=figure('Name','Confusion matrix', 'MenuBar', 'none', 'NumberTitle', 'off');
    if (~isempty(test_path))&(size(test_data,2)>cov_dim),
        ztest=fhme(bnet2, nodes_info, test_d, size(test_d,1));
        [Htest, testRate]=confmat(ztest, test_t);   % CM on the test set
        subplot(1,2,1);
    end
    plotmat(Htrain,'b','k',12)
    tick=[0.5:1:(0.5+nodes_info(2,end)-1)];
    set(gca,'XTick',tick)
    set(gca,'YTick',tick)
    grid('off')
    ylabel('True')
    xlabel('Prediction')
    title(['Confusion Matrix: training set (' num2str(trainRate(1)) '%)'])
    if (~isempty(test_path))&(size(test_data,2)>cov_dim),
        subplot(1,2,2)
        plotmat(Htest,'b','k',12)
        set(gca,'XTick',tick)
        set(gca,'YTick',tick)
        grid('off')
        ylabel('True')
        xlabel('Prediction')
        title(['Confusion Matrix: test set (' num2str(testRate(1)) '%)'])
    end
    disp(' ')
    disp('(Press any key to continue)');
    pause
end
%-----------------------------------------------------------------------------------
% Regression & Classification problem: calculate the predictions & plot the LL trace
%-----------------------------------------------------------------------------------
train_result=fhme(bnet2,nodes_info,train_d,size(train_d,1));
if ~isempty(test_path),
    test_result=fhme(bnet2,nodes_info,test_d,size(test_d,1));
end
fh3=figure('Name','Log-likelihood trace', 'MenuBar', 'none', 'NumberTitle', 'off')
plot(LL2,'-ro',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[1 1 0],...
                'MarkerSize',4)
title('Log-likelihood trace')
%-----------------------------------------------------------------------------------
% Regression & Classification problem: save the predictions
%-----------------------------------------------------------------------------------
clc
disp('------------------------------------------------------------------');
disp('                           Save the results                       ');
disp('------------------------------------------------------------------');
disp(' ')
%-----------------------------------------------------------------------------------
save_quest_m=input('Do you want to save the HME model (Y/N)? [Y default]: ', 's');
if isempty(save_quest_m),
    save_quest_m='Y';
end
if ~findstr(save_quest_m, ['Y', 'N']), error('Invalid input'); end
if save_quest_m=='Y',
    disp(' ');
    m_save=input('Insert the complete path for save the HME model (.mat):\n >> ', 's');
    if isempty(m_save), error('You must specify a path!'); end
    save(m_save, 'bnet2');
end
%-----------------------------------------------------------------------------------    
disp(' ')
save_quest=input('Do you want to save the HME predictions (Y/N)? [Y default]: ', 's');
disp(' ')
if isempty(save_quest),
    save_quest='Y';
end
if ~findstr(save_quest, ['Y', 'N']), error('Invalid input'); end
if save_quest=='Y',
    tr_save=input('Insert the complete path for save the training data prediction (.mat):\n >> ', 's');    
    if isempty(tr_save), error('You must specify a path!'); end
    save(tr_save, 'train_result');  
    if ~isempty(test_path),
        disp(' ')
        te_save=input('Insert the complete path for save the test data prediction (.mat):\n >> ', 's');
        if isempty(te_save), error('You must specify a path!'); end
        save(te_save, 'test_result');
    end
end
clc
disp('----------------------------------------------------');
disp('                      B  Y  E !                     ');
disp('----------------------------------------------------');
pause(2)
%clear 
clc
