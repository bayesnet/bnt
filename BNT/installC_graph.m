% These C functions were written by Ilya Shpitser.

if isunix
  mex -c elim.c;
  mex -c cell.c;
  mex -c map.c;
  mex -DUNIX  best_first_elim_order.c  elim.o  cell.o  map.o;
  mex -DUNIX  triangulate.c  elim.o  cell.o  map.o;
else
  mex -c elim.c;
  mex -c cell.c;
  mex -c map.c;
  mex  best_first_elim_order.c  elim.obj  cell.obj  map.obj;
  mex  triangulate.c  elim.obj  cell.obj  map.obj;
end

