function S = export_dnet(bnet, file, proba)
% filepath = export_dnet(bnet, 'filename', includeparameters)
%  filename and includeparameters ([0] or 1) are optional
%
% Exports BNT bnets to Netica dnet
%
% Written by Francois.Olivier.C.H@gmail.com
%
% Informations could be found here http://www.norsys.com/downloads/
%
% Version 80806
% Supports only discret bayesian network

if nargin==1, file=['dnet' datestr(now,'-yymmdd-HHMMSS')]; end
if nargin<3, proba=0; end

% generating filename
if length(file)>4,
  if prod(double(file((end-4):end)~='.dnet')), file=[file '.dnet']; end
  name = file(1:end-5);
else
  name = file;
  file = [file '.dnet'];
end

% generating node names 1:N if non existant
N=length(bnet.dag);
if isempty(bnet.names),
  for i=1:N,
    keys{i}=['' num2str(i) ''];
    vals{i}=i;
  end
  bnet.names = assocarray(keys, vals);
end

% header of the file
fid = fopen(file, 'w');
fprintf(fid, '// ~->[DNET-1]->~\n\n');
fprintf(fid, '// exported from the Bayes Net Toolbox with export_dnet function \n');
fprintf(fid, '// please report bugs to francois.olivier.c.h@gmail.com\n');
if proba, fprintf(fid, '// Take care !  Parents'' order isn''t the same in probability table comments\n'); end
fprintf(fid, ['\nbnet ' name ' {']);

% main loop
for node = 1:N
  name = get_key(bnet.names,node);

  fprintf(fid, ['\nnode ' name ' {']);

  fprintf(fid, '\n\tkind = NATURE;');

  fprintf(fid, '\n\tdiscrete = TRUE;');

% states names are x1,x2:xsize(node)
  fprintf(fid, '\n\tstates = (');
  for l=1:bnet.node_sizes(node),
    fprintf(fid, 'x%d',l);
    if l~=bnet.node_sizes(node), fprintf(fid, ', '); end
  end, fprintf(fid, ');');

% declare parent in counter order to be coherent with prob section
  fprintf(fid, '\n\tparents = (');
  par = bnet.parents{node};
  for l=length(par):-1:1
    fprintf(fid, '%s',get_key(bnet.names,par(l)));
    if l~=1, fprintf(fid, ', '); end
  end, fprintf(fid, ');');

% fill probs if proba==1
    if proba,
        fprintf(fid, '\n\tprobs =\n');
        % inits
        fprintf(fid, '\t//\t');
        for l = 1:bnet.node_sizes(node),
            fprintf(fid, 'x%d\t', l);
        end
            fprintf(fid, '\t//');
        if ~isempty(par),
            for l = unique(par), %right order this time because of the way both BNT and netica work
                fprintf(fid, '\t%s', get_key(bnet.names,l));
            end
        end

        % opens tab
        fprintf(fid, '\n\t');
        for l = unique(par),
            fprintf(fid, '(');
        end

        % fullfils probs
        CPT = CPT_from_bnet(bnet);
        CPT = CPT{node};
        CPT=CPT(:);            % good order whatever the node size ????
        if isempty(par),
            fprintf(fid, '\t');
            for i=1:length(CPT)-1,
                fprintf(fid, '%1.4f, ', CPT(i));
            end
            fprintf(fid, '%1.4f);',CPT(end));
        else                   % if there are parents
            endi=0;
            parsiz = prod(bnet.node_sizes([par]));
            parentstates = ones(1,length(par)); parentstates(end)=0;
            for i=1:parsiz
                % prints probas
                fprintf(fid, '(\t');
                for j = 1:bnet.node_sizes(node)
                    prob = CPT(i+parsiz*(j-1));
                    if j~=bnet.node_sizes(node), fprintf(fid, '%1.4f, ',prob);
                    else fprintf(fid, '%1.4f',prob); end
                end

                % closes parenthesis if needed
                endii=endi;
                if i~=parsiz,
                    while endi>0, fprintf(fid, ')');endi=endi-1; end
                    fprintf(fid, '),\t//\t');
                else % close tab
                    for l = unique(par), fprintf(fid, ')'); end
                    fprintf(fid, ');\t//\t');
                end

                % prints node states
                res=i;
                for l=1:length(par)-1
                    if mod(i+1, prod(bnet.node_sizes(par(1:l))))==0, endi=endi+1; end  % counts parenthesis

                    resaff = mod(res, bnet.node_sizes(l));
                    if resaff==0, resaff=bnet.node_sizes(par(l)); end
                    res = div(res-1, bnet.node_sizes(par(l)))+1;

                    fprintf(fid, '%d\t',resaff);
                end

                resaff = div(i-1, prod(bnet.node_sizes(par(1:end-1))))+1;
                fprintf(fid, '%d',resaff);

                % opens parenthesis if needed
                if i~=parsiz, fprintf(fid, '\n'); end
                fprintf(fid, '\t');
                if i~=parsiz, while endii>0, fprintf(fid, '(');endii=endii-1; end, end
            end
        end
    end
    fprintf(fid, '\n\t};');
end

% closes file
fprintf(fid,'\n};\n');
fclose(fid);

% outputs string
S = [pwd '/' file];
