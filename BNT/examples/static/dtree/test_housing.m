% Here the training data is adapted from UCI ML repository, 'housing' data
% Input variables: 12 continous, one binary
% Ouput variables: continous
% The testing result trace is in the end of this script, it is same to the graph in page 219 of 
% Leo Brieman etc. 1984 book titled "Classification and regression trees".

dtreeCPD=tree_CPD;

% load data
fname = fullfile(BNT_HOME, 'examples', 'static', 'uci_data', 'housing', 'housing.data');
data=load(fname);
data=data';
data=transform_data_into_bnt_format(data,[1:3,5:14]); 

% learn decision tree from data 
ns=1*ones(1,14);
ns(4)=2;
dtreeCPD1=learn_params(dtreeCPD,1:14,data,ns,[1:3,5:14],'stop_cases',5,'min_gain',0.006); 

% evaluate on data
[score,outputs]=evaluate_tree_performance(dtreeCPD1,1:14,data,ns,[1:3,5:14]);
fprintf('Mean square deviation (using regression tree to predict) in old training data %6.3f\n',score);


% show decision tree using graphpad
% It should be easy, but still not implemented



% >> test_housing
% Create node 1 split at 6 gain  38.2205 Th 6.939000e+000. Mean  22.5328 Cases 506
% Create node 2 split at 13 gain  14.4503 Th 1.437000e+001. Mean  19.9337 Cases 430
% Create node 3 split at 8 gain   4.9809 Th 1.358000e+000. Mean  23.3498 Cases 255
% Create node 4 split at 1 gain   0.7722 Th 1.023300e+001. Mean  45.5800 Cases 5
% Create leaf node(samevalue) 5. Mean  50.0000 Std   0.0000 Cases 4 
% Add subtree node 5 to 4. #nodes 5
% Create leaf node(samevalue) 6. Mean  27.9000 Std   0.0000 Cases 1 
% Add subtree node 6 to 4. #nodes 6
% Add subtree node 4 to 3. #nodes 6
% Create node 7 split at 6 gain   2.8497 Th 6.540000e+000. Mean  22.9052 Cases 250
% Create node 8 split at 13 gain   0.5970 Th 7.560000e+000. Mean  21.6297 Cases 195
% Create leaf node(nogain) 9. Mean  23.9698 Std   1.7568 Cases 43 
% Add subtree node 9 to 8. #nodes 9
% Create leaf node(nogain) 10. Mean  20.9678 Std   2.8242 Cases 152 
% Add subtree node 10 to 8. #nodes 10
% Add subtree node 8 to 7. #nodes 10
% Create leaf node(nogain) 11. Mean  27.4273 Std   3.4512 Cases 55 
% Add subtree node 11 to 7. #nodes 11
% Add subtree node 7 to 3. #nodes 11
% Add subtree node 3 to 2. #nodes 11
% Create node 12 split at 1 gain   2.2467 Th 6.962150e+000. Mean  14.9560 Cases 175
% Create node 13 split at 5 gain   0.5172 Th 5.240000e-001. Mean  17.1376 Cases 101
% Create leaf node(nogain) 14. Mean  20.0208 Std   3.0672 Cases 24 
% Add subtree node 14 to 13. #nodes 14
% Create leaf node(nogain) 15. Mean  16.2390 Std   2.9746 Cases 77 
% Add subtree node 15 to 13. #nodes 15
% Add subtree node 13 to 12. #nodes 15
% Create node 16 split at 5 gain   0.6133 Th 6.050000e-001. Mean  11.9784 Cases 74
% Create leaf node(nogain) 17. Mean  16.6333 Std   4.5052 Cases 12 
% Add subtree node 17 to 16. #nodes 17
% Create leaf node(nogain) 18. Mean  11.0774 Std   3.0090 Cases 62 
% Add subtree node 18 to 16. #nodes 18
% Add subtree node 16 to 12. #nodes 18
% Add subtree node 12 to 2. #nodes 18
% Add subtree node 2 to 1. #nodes 18
% Create node 19 split at 6 gain   6.0493 Th 7.420000e+000. Mean  37.2382 Cases 76
% Create node 20 split at 1 gain   1.9900 Th 7.367110e+000. Mean  32.1130 Cases 46
% Create node 21 split at 8 gain   0.6273 Th 1.877300e+000. Mean  33.3488 Cases 43
% Create leaf node(samevalue) 22. Mean  45.6500 Std   6.1518 Cases 2 
% Add subtree node 22 to 21. #nodes 22
% Create leaf node(nogain) 23. Mean  32.7488 Std   3.5690 Cases 41 
% Add subtree node 23 to 21. #nodes 23
% Add subtree node 21 to 20. #nodes 23
% Create leaf node(samevalue) 24. Mean  14.4000 Std   3.7363 Cases 3 
% Add subtree node 24 to 20. #nodes 24
% Add subtree node 20 to 19. #nodes 24
% Create node 25 split at 1 gain   1.1001 Th 2.733970e+000. Mean  45.0967 Cases 30
% Create leaf node(nogain) 26. Mean  45.8966 Std   4.4005 Cases 29 
% Add subtree node 26 to 25. #nodes 26
% Create leaf node(samevalue) 27. Mean  21.9000 Std   0.0000 Cases 1 
% Add subtree node 27 to 25. #nodes 27
% Add subtree node 25 to 19. #nodes 27
% Add subtree node 19 to 1. #nodes 27
% Mean square deviation (using regression tree to predict) in old training data  9.405
% 


