function pretty_print_hhmm_parse(mpe, Qnodes, Fnodes, Onode, alphabet)
% function pretty_print_hhmm_parse(mpe, Qnodes, Fnodes, Onode, alphabet)
%
% mpe(i,t) is the most probable value of node i at time t
% Qnodes(1:D), Fnodes = [F2 .. FD], Onode contain the node ids
% alphabet(i) is the i'th output symbol, or [] if don't want displayed

T = size(mpe,2);
ncols = 20;
t1 = 1; t2 = min(T, t1+ncols-1);
while (t1 < T)
  %fprintf('%d:%d\n', t1, t2);
  if iscell(mpe)
    print_block_cell(mpe(:,t1:t2), Qnodes, Fnodes, Onode, alphabet, t1);
  else
    print_block(mpe(:,t1:t2), Qnodes, Fnodes, Onode, alphabet, t1);
  end
  fprintf('\n\n');
  t1 = t2+1; t2 = min(T, t1+ncols-1);
end

%%%%%%

function print_block_cell(mpe, Qnodes, Fnodes, Onode, alphabet, start)

D = length(Qnodes);
T = size(mpe, 2);
fprintf('%3d ', start:start+T-1); fprintf('\n');
for d=1:D
  for t=1:T
    if (d > 1) & (mpe{Fnodes(d-1),t} == 2)
      fprintf('%3d|', mpe{Qnodes(d), t});
    else
      fprintf('%3d ', mpe{Qnodes(d), t});
    end
  end
  fprintf('\n');
end
if ~isempty(alphabet)
  a = cell2num(mpe(Onode,:));
  %fprintf('%3c ', alphabet(mpe{Onode,:}));
  fprintf('%3c ', alphabet(a))
  fprintf('\n');
end


%%%%%%

function print_block(mpe, Qnodes, Fnodes, Onode, alphabet, start)

D = length(Qnodes);
T = size(mpe, 2);
fprintf('%3d ', start:start+T-1); fprintf('\n');
for d=1:D
  for t=1:T
    if (d > 1) & (mpe(Fnodes(d-1),t) == 2)
      fprintf('%3d|', mpe(Qnodes(d), t));
    else
      fprintf('%3d ', mpe(Qnodes(d), t));
    end
  end
  fprintf('\n');
end
if ~isempty(alphabet)
  fprintf('%3c ', alphabet(mpe(Onode,:)));
  fprintf('\n');
end
