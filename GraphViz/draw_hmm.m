function draw_hmm(A, varargin)
% DRAW_HMM Make a picture of the HMM using dotty
% function draw_hmm(A, ...)
%
% For details on dotty, see http://www.research.att.com/sw/tools/graphviz
%
% If A(i,j) > thresh, we draw and arc from state i to state j.
%
% Optional arguments (name/value pairs) [default]
%
% thresh - [1e-1]
% obsprob - If B(i,o) > 0, we include "o" in the name of state i.
%     e.g., if state 5 emits 1,3,7, its label becomes "5: 1 3 7".
% startprob - ifstartprob(i) > 0, the state name will be prefixed with "+".
% endprob - if endprob(i) > 0, the state name will be appended with "-".
% filename - if [], we write to 'tmp.dot', convert this to 'tmp.ps'
%  using 'dot -Tps tmp.dot -o tmp.ps', and then call ghostview  to display the result.
%  dot and gv must be on your system path.
%  If filename ~= [], we just generate the dot file, and do not
%  convert it to postscript or call ghostview.

[thresh, B, startprob, endprob, filename] = ...
    process_options(varargin, 'thresh', 1e-1, 'obsprob', [], 'startprob', [], 'endprob', [], ...
		    'filename', []);

Q = length(A);

arclabel = cell(Q,Q);
G = zeros(Q,Q);
for i=1:Q
  for j=1:Q
    if A(i,j) < thresh
      arclabel{i,j} = '';
    else
      G(i,j) = 1;
      arclabel{i,j} = sprintf('%5.3f', A(i,j));
    end
  end
end


nodelabel = cell(1,Q);
for i=1:Q
  % annotate start/stop states
  if ~isempty(startprob) & ~approxeq(startprob(i), 0)
    start = '+';
  else
    start = '';
  end
  if ~isempty(endprob) & ~approxeq(hmm.endprob(i), 0)
    stop = '-';
  else
    stop = '';
  end
  label = sprintf('%s%d%s :', start, i, stop);

  if ~isempty(B)
    output_label = mk_output_label(B);
    label = strcat(label, output_label);
  end
  
  nodelabel{i} = label;
end


if isempty(filename)
  filename = 'tmp.dot';
  %mkdot(G, filename, arclabel, nodelabel)  
  graph_to_dot(G, 'filename', filename, 'arc_label', arclabel, 'node_label', nodelabel);
  fprintf('converting from .ps to .dot\n')
  !dot -Tps tmp.dot -o tmp.ps
  !gv tmp.ps &
else
  graph_to_dot(G, 'filename', filename, 'arc_label', arclabel, 'node_label', nodelabel);
  %mkdot(G, filename, arclabel, nodelabel)  
end


%%%%%%%%%

function label = mk_output_label(B)

[Q O] = size(B);
label = '';

if 0
  % print most probable symbols
  for i=1:Q
    m = max(B(i,:));
    ndx = find(abs(B(i,:) - repmat(m,1,O)) < 1e-2);
    %ndx = find(B(i,:)==m);
    %label = sprintf('%d,', ndx);
  end
end

if 0
  % print prob distrib over all symbols 
  for o=1:O
    if approxeq(B(i,o), 0)
      %
    else
      label = strcat(label, sprintf('%d(%3.2f),', o, B(i,o)));
    end
  end
end

if 1
  % print all non-zero symbols
  chars = ['a' 'b' 'c'];
  for o=1:O
    if approxeq(B(i,o), 0)
      %
    else
      label = strcat(label, sprintf('%s', chars(o)));
    end
  end
end
