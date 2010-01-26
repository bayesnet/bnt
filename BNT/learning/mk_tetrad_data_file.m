function mk_tetrad_data_file(filename, samples, sig)
% MK_TETRAD_DATA_FILE Make a file containing raw discrete data for input to TETRAD
% mk_tetrad_data_file(filename, samples, sig)
%
% samples(i,j) is the value for case i, variable j
% The resulting file can be used for the 'build' part of Tetrad.
% For details on tetrad, see hss.cmu.edu/html/departments/philosophy/TETRAD/tetrad.html

[nsamples N] = size(samples);

fid = fopen(filename, 'w');
fprintf(fid, '/Raw\n');
fprintf(fid, '%d\n', nsamples);
for i=1:N
  fprintf(fid, 'x%d ', i);
end
fprintf(fid, '\n');
for i=1:nsamples
  fprintf(fid, '%d ', samples(i,:)-1); % tetrad counts from 0
  fprintf(fid, '\n');
end
%fprintf(fid, '/Knowledge\n');
%fprintf(fid, 'Significance %4.2f\n', sig);
fclose(fid);

