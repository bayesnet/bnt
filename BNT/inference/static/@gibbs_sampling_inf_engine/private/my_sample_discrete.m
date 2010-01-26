function M = my_sample_discrete(prob)
% A faster version that calls a c subfunction.  Will update one
% day to have r and c parameters as well

R = rand (1,1);
M = sample_single_discrete(R, prob);

