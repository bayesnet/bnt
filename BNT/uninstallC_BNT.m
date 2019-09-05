dirs = {'potentials/Tables', ...
	'CPDs/@discrete_CPD', ...
      'inference/static/@jtree_sparse_inf_engine', ...
      'inference/static/@gibbs_sampling_inf_engine/private'};

BNT_HOME = '/home/ai2/murphyk/matlab/FullBNT'; % edit this
%global BNT_HOME

for d=1:length(dirs)
  f = fullfile(BNT_HOME, 'BNT', dirs{d});
  fprintf('removing Cmex files from %s\n', f);
  cd(f)
  delete *.mex*
  delete *.dll
  delete *.obj
  delete *.o
end
