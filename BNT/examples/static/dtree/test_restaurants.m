% Here the training data is adapted from Russell95 book. See restaurant.names for description.
% (1) Use infomation-gain as the split testing score, we get the the same decision tree as the book Russell 95 (page 537),
% and the Gain(Patrons) is 0.5409, equal to the result in Page 541 of Russell 95. (see below output trace)
% (Note: the dtree in that book has small compilation error, the Type node is from YES of Hungry node, not NO.)
% (2) Use gain-ratio (Quilan 93), the splitting defavorite attribute with more values. (e.g. the Type attribute here)

dtreeCPD=tree_CPD;

% load data
fname = fullfile(BNT_HOME, 'examples', 'static', 'uci_data', 'restaurant', 'restaurant.data');
data=load(fname);
data=data';

%make the data be BNT compliant (values for discrete nodes are from 1-n, here n is the node size)
  % e.g. if the values are [0 1 6], they must be mapping to [1 2 3]
%data=transform_data(data,'tmp.dat',[]); %here no cts nodes

% learn decision tree from data 
ns=2*ones(1,11);
ns(5:6)=3;
ns(9:10)=4;
dtreeCPD1=learn_params(dtreeCPD,1:11,data,ns,[]);

% evaluate on data
[score,outputs]=evaluate_tree_performance(dtreeCPD1,1:11,data,ns,[]);
fprintf('Accuracy in training data %6.3f\n',score);

% show decision tree using graphpad



% --------------------------Output trace: using Information-Gain------------------------------
% The splits are Patron, Hungry, Type, Fri/Sat
% *********************************
% Create node 1 split at 5 gain 0.5409 Th 0. Class 1 Cases 12 Error 6 
% Create leaf node(onecla) 2. Class 1 Cases 2 Error 0 
% Add subtree node 2 to 1. #nodes 2
% Create leaf node(onecla) 3. Class 2 Cases 4 Error 0 
% Add subtree node 3 to 1. #nodes 3
% Create node 4 split at 4 gain 0.2516 Th 0. Class 1 Cases 6 Error 2 
% Create leaf node(onecla) 5. Class 1 Cases 2 Error 0 
% Add subtree node 5 to 4. #nodes 5
% Create node 6 split at 9 gain 0.5000 Th 0. Class 1 Cases 4 Error 2 
% Create leaf node(nullset) 7. Father 6 Class 1
% Create node 8 split at 3 gain 1.0000 Th 0. Class 1 Cases 2 Error 1 
% Create leaf node(onecla) 9. Class 1 Cases 1 Error 0 
% Add subtree node 9 to 8. #nodes 9
% Create leaf node(onecla) 10. Class 2 Cases 1 Error 0 
% Add subtree node 10 to 8. #nodes 10
% Add subtree node 8 to 6. #nodes 10
% Create leaf node(onecla) 11. Class 2 Cases 1 Error 0 
% Add subtree node 11 to 6. #nodes 11
% Create leaf node(onecla) 12. Class 1 Cases 1 Error 0 
% Add subtree node 12 to 6. #nodes 12
% Add subtree node 6 to 4. #nodes 12
% Add subtree node 4 to 1. #nodes 12
% ********************************
% 
% Note:
% ***Create node 4 split at 4 gain 0.2516 Th 0. Class 1 Cases 6 Error 2 
% This mean we create a new node number 4, it is splitting at the attribute 4, and info-gain is 0.2516, 
% "Th 0" means threshhold for splitting continous attribute, "Class 1" means the majority class at node 4 is 1,
% and "Cases 6" means it has 6 cases attached to it, "Error 2" means it has two errors if changing the class lable of 
% all the cases in it to the majority class.
% *** Add subtree node 12 to 6. #nodes 12
% It means we add the child node 12 to node 6.
% *** Create leaf node(onecla) 10. Class 2 Cases 1 Error 0 
% here 'onecla' means all cases in this node belong to one class, so no need to split further. 
%      'nullset' means no training cases belong to this node, we use its parent node majority class as its class
% 
% 
% 
% ---------------Output trace: using GainRatio-----------------------
% The splits are Patron, Hungry, Fri/Sat, Price
% 
% 
% Create node 1 split at 5 gain 0.3707 Th 0. Class 1 Cases 12 Error 6 
% Create leaf node(onecla) 2. Class 1 Cases 2 Error 0 
% Add subtree node 2 to 1. #nodes 2
% Create leaf node(onecla) 3. Class 2 Cases 4 Error 0 
% Add subtree node 3 to 1. #nodes 3
% Create node 4 split at 4 gain 0.2740 Th 0. Class 1 Cases 6 Error 2 
% Create leaf node(onecla) 5. Class 1 Cases 2 Error 0 
% Add subtree node 5 to 4. #nodes 5
% Create node 6 split at 3 gain 0.3837 Th 0. Class 1 Cases 4 Error 2 
% Create leaf node(onecla) 7. Class 1 Cases 1 Error 0 
% Add subtree node 7 to 6. #nodes 7
% Create node 8 split at 6 gain 1.0000 Th 0. Class 2 Cases 3 Error 1 
% Create leaf node(onecla) 9. Class 2 Cases 2 Error 0 
% Add subtree node 9 to 8. #nodes 9
% Create leaf node(nullset) 10. Father 8 Class 2
% Create leaf node(onecla) 11. Class 1 Cases 1 Error 0 
% Add subtree node 11 to 8. #nodes 11
% Add subtree node 8 to 6. #nodes 11
% Add subtree node 6 to 4. #nodes 11
% Add subtree node 4 to 1. #nodes 11
% 
% 
