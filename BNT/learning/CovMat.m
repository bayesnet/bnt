function [CovMatrix, obs, varfields] = CovMat(filename,row_cols)
%[CovMatrix, obs, varfields] = CovMat(filename,row_cols)
%% generates a Covariance Matrix from a file of data consisting of N columns of M data rows
%%      filename        string name (with path and extension) of file to open
%%      row_cols        Number_of_converstions_per_row  (turns into  [3 inf])
%% Return
%%      CovMatrix       Covariance matrix
%%      obs             Number of observations read in
%%      varfields       Labels of the variables see filename structure below
%%
%%  Filename structure: 
%%      Comma separated, starting with the variable labels, then the data in rows.
%%    filename test.txt consists of:
%%
%%      Earthquake,Burglar,Radio,Alarm,Call
%%      1,2,3,4,5
%%      11,22,33,44,55
%%      . . .
%%
%% Example call: 
%%  [cvmat numdat lables] = CovMat('test.txt',5);
%%
%%      Returns Covariance matrix, number of date rows and variable field names
%% Gary R. Bradski   7/2002

fmtstr = '%f';
for i = 2:row_cols
    fmtstr = strcat(fmtstr,',%f');
end

%% load data
fidCov = fopen(filename,'r');

varfields = fgetl(fidCov);
Corx = fscanf(fidCov,fmtstr,[row_cols inf]);
Corx= Corx';
[obs bla] = size(Corx);
CovMatrix = cov(Corx); 
fclose(fidCov);