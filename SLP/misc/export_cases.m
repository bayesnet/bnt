function S = export_cases(data, names, file, misv)
% filepath = export_cases(data, names, 'filename', missing_value)
%  filename and missing_value [default [] if iscell(data) or -9999 if not] are optional
% 
% Exports BNT datasets to Netica cases (* for missing data)
%
% Written by Francois.Olivier.C.H@gmail.com
% 
% Informations could be found here http://www.norsys.com/downloads/
%
% version 80807

% inits
if nargin==1, file=['dnet' datestr(now,'-yymmdd-HHMMSS')]; end
if nargin<2, error('Variable Names needed'); end
if length(file)>5,
  if prod(double(file((end-3):end)~='.cas')), file=[file '.cas']; end
  name = file(1:end-5);
else
  name = file;
  file = [file '.cas'];
end
if nargin<4, misv=-9999; end
if iscell(data), data = bnt_to_mat(data,misv); end

[N m]=size(data);
if length(names)~=N, error('Sizes must be the same'); end

% header of the file
fid = fopen(file, 'w');
fprintf(fid, '// exported from the Bayes Net Toolbox with export_cases function \n');
fprintf(fid, '// please report bugs to francois.olivier.c.h@gmail.com\n\n');

% write names
for i=1:N
  fprintf(fid, '%s\t',names{i});
end

% exports
for l=1:m,
  fprintf(fid, '\n');
  for i=1:N,
    if data(i,l)~=misv,
      fprintf(fid, '%s',['x' num2str(data(i,l))]);
    else
      fprintf(fid, '%s','*');
    end
    if i<N, fprintf(fid, '\t'); end
  end
end

% closes file
fprintf(fid,'\n');
fclose(fid);

% outputs string
S = [pwd '/' file];
