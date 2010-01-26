function seq = sample_dbn(bnet, varargin)
% SAMPLE_DBN Generate a random sequence from a DBN.
% seq = sample_dbn(bnet, ...)
%
% seq{i,t} contains the values of the i'th node in the t'th slice.
%
% Optional arguments:
%
% length - length of sequence  to be generated (can also just use sample_dbn(bnet,T))
% stop_test - name of a function which is used to decide when to stop;
%   This will be called as feval(stop_test, seq(:,t))
%   i.e., stop_test is passed a cell array containing all the nodes in the current slice.   
% evidence - initial evidence; if evidence{i,t} is non-empty, this node won't be sampled.

args = varargin;
nargs = length(args);

if (nargs == 1) & ~isstr(args{1})
  % Old syntax: sample_dbn(bnet, T)
  T = args{1};
else
  % get length
  T = 1;
  for i=1:2:nargs
    switch args{i},
     case 'length',      T = args{i+1}; 
     case 'evidence',    T = size(args{i+1}, 2);
    end
  end
end

ss = length(bnet.intra);
% set default arguments
seq = cell(ss, T);
stop_test = [];
for i=1:2:nargs
  switch args{i},
   case 'evidence',    seq = args{i+1}; % initialise observed nodes
   case 'stop_test',  stop_test = args{i+1};
  end
end

t = 1;
for i=1:ss
  if ~isempty(stop_test) | isempty(seq{i,t})
    ps = parents(bnet.dag, i);
    e = bnet.equiv_class(i,1);
    pvals = seq(ps);
    seq{i,t} = sample_node(bnet.CPD{e}, pvals);
    %fprintf('sample i=%d,t=%d,val=%d,ps\n', i, t, seq(i,t)); pvals(:)'
  end
end
t = 2;
done = 0;
while ~done
  for i=1:ss
    if ~isempty(stop_test) | isempty(seq{i,t})
      ps = parents(bnet.dag, i+ss) + (t-2)*ss;
      e = bnet.equiv_class(i,2);
      pvals = seq(ps);
      seq{i,t} = sample_node(bnet.CPD{e}, pvals);
      %fprintf('sample i=%d,t=%d,val=%d,ps\n', i, t, seq(i,t)); pvals(:)'
    end
  end
  if ~isempty(stop_test)
    done = feval(stop_test, seq(:,t));
  else
    if t==T
      done = 1;
    end
  end
  t = t + 1;
end
