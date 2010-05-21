function [data, comp_data, bnet_miss, taux, bnet_orig, notok, d] = gener_data_from_bnet_miss(bnet_miss, m, base_proba ,v, testdata)
% [data, comp_data, bnet_miss, taux, bnet_orig, notok] = gener_data_from_bnet_miss(bnet_miss, m, base_proba ,v, aretestdata)
%
%   bnet_miss : see gener_[MCAR or MAR]_net function
%   m : the length of the dataset
%   if base_proba==0 or does not exist the x2_test will be passed
%   v==1 to enter the verbose mode
%   aretestdata==1 to always build the same dataset <-- rand('state',0)
%
% Francois.Olivier.C.H@gmail.com

% Initialisation

if nargin<5, testdata = 0; end
if nargin<4, v = 0; end

N2 = length(bnet_miss.dag); 
if mod(N2,3)~=0, error('The number of nodes must be even in bnet_miss'); end
N = length(bnet_miss.dag)/3;

if nargin<3, base_proba=0; end
if nargin<2, error('Not enougth parameters'); end

% CHOOSE THE TEST POWER (only affect 'notok' value)
                    chi2_0_1_1fd  =  2.705 ;
                   chi2_0_05_1fd  =  3.841 ;
                   chi2_0_01_1fd  =  6.635 ;
                   chi2_0_001_1fd = 10.827 ;
                  chi2_0_0001_1fd = 15.137 ;
                                            choice = chi2_0_001_1fd;
notok=0;
clear chi2_0_1_1fd chi2_0_05_1fd chi2_0_01_1fd chi2_0_001_1fd chi2_0_0001_1fd

% Recovering bnet_orig

dag = bnet_miss.dag(1:N,1:N);
bnet_orig = mk_bnet(dag, bnet_miss.node_sizes(1:N));
CPT = CPT_from_bnet(bnet_miss, 0);
for i=1:N, bnet_orig.CPD{i} = tabular_CPD(bnet_orig, i, CPT{i}); end

% Generation of complete data

if testdata, rand('state',0); randn('state',0); end

data = cell(N,m);
for l = 1:m, data(:,l) = sample_bnet(bnet_orig); end
fprintf('Complete data have been created.');

% Generation of missing array

miss_array = cell(3*N,m);
vide = cell(1,2*N); l= 1;
while l <= m, 
  ev(1:N) = data(:,l); ev(N+1:3*N) = vide;
  miss_array(:,l) = sample_bnet(bnet_miss, 'evidence', ev); 
  % apply simple rule of missingness
  ev2 = cell2mat(miss_array(N+1:2*N, l));
  % verification that we have not a completly missing sample
  ev2 = 3-ev2; 
  if prod(ev2)==1, 
    if v, fprintf(' - %d, one completly missing sample removed', l);end
  else l=l+1;
  end 
  if v, if mod(l,250)==0, fprintf('\n - %d',l); end,  end
end
fprintf('\n');

% Generation of incomplete dataset

       %%   TO REPLACE THE CELL ARRAY FOR ouput data
       %%   WITH A MATRIX WITH A SPECIAL CASE (size+1)
       %%   FOR MISSING DATA, SIMPLY REPLACE 1 by 0
if 1,  %%   HERE
  miss_array = bnt_to_mat(miss_array(N+1:2*N, :)); 
  miss_array = 2-miss_array; 
  data = bnt_to_mat(data);
  comp_data = data; 
  data = data.*miss_array; 
  data = mat_to_bnt(data, 0); 
else
  comp_data = bnt_to_mat(miss_array(1:N,:),0);
  data = bnt_to_mat(miss_array(2*N+1:3*N,:),0);
  miss_array = bnt_to_mat(miss_array(N+1:2*N, :)); 
  miss_array = 2-miss_array; 
end
fprintf('Incomplete dataset have been created.\n');

% Verification of the Rate of missing data

if base_proba,
    [XX, YY]=find(miss_array==0);
    nbr_miss = length(YY);
    taux = nbr_miss/N/m;
    if v, fprintf('There is %2.2f percent of missing data\n', round(taux*10000)/100); end

    % Khi2 test between taux and base_proba for m*N
    toto = m*N;
    d = ((nbr_miss-base_proba*toto)^2)/(base_proba*toto) + (((toto-nbr_miss)-(1-base_proba)*toto)^2)/((1-base_proba)*toto);
    if d>choice, 
    fprintf('THE DATASET DO NOT RESPECT %2.1f%% OF MISSING DATA (%2.1f%%, Khi2 : %2.1f > %2.1f)\n', round(base_proba*10000)/100, round(taux*10000)/100, d, choice); 
    notok = 1;
    end
end
