function mk_collage_from_clqs(dir, cliques)

% For use with mk_ps_from_clqs.
% This generates a latex file that glues all the .ps files
% into one big figure.

cd(dir)
C = length(cliques);

ncols = 4;
width = 1.5;
fid = fopen('collage.tex', 'w');
fprintf(fid, '\\documentclass{article}\n');
fprintf(fid, '\\usepackage{psfig}\n');
fprintf(fid, '\\begin{document}\n');
fprintf(fid, '\\centerline{\n');
fprintf(fid, '\\begin{tabular}{');
for col=1:ncols,  fprintf(fid, 'c'); end
fprintf(fid, '}\n');
c = 1;
for row = 1:floor(C/ncols)
  for col=1:ncols-1
    fname = sprintf('%s/clq%d.ps', dir, c);
    fprintf(fid, '\\psfig{file=%s,width=%3fin} & \n', fname, width);
    c = c + 1;
  end
  fname = sprintf('%s/clq%d.ps', dir, c);
  fprintf(fid, '\\psfig{file=%s,width=%3fin} \\\\ \n', fname, width);
  c = c + 1;
end
% last row
while (c <= C)
  fname = sprintf('%s/clq%d.ps', dir, c);
  fprintf(fid, '\\psfig{file=%s,width=%3fin} & \n', fname, width);
  c = c + 1;
end
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '}\n');
fprintf(fid, '\\end{document}');
fclose(fid);

!latex collage.tex &
!dvips -o collage.ps collage.dvi &
!ghostview collage.ps &
