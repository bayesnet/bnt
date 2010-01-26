function [intra, inter, Qnodes, Fnodes, Onode] = mk_hhmm_topo(D, all_Q_to_Qs, Ops, F1)
% MK_HHMM_TOPO Make Hierarchical HMM topology
% function [intra, inter, Qnodes, Fnodes, Onode] = mk_hhmm_topo(D, all_Q_to_Qs, Ops, F1)
%
% D is the depth of the hierarchy
% If all_Q_to_Qs = 1, level i connects to all levels below, else just to i+1 [0]
% Ops are the Q parents of the observed node [Qnodes(end)]
% If F1=1, level 1 can finish (restart), else there is no F1->Q1 arc [0]

Qnodes = 1:D;

if nargin < 2, all_Q_to_Qs = 1; end
if nargin < 3, Ops = Qnodes(D); end
if nargin < 4, F1 = 0; end

if F1
  Fnodes = 2*D:-1:D+1; % must number from bottom to top
  Onode = 2*D+1;
  ss = 2*D+1;
else
  Fnodes = [-1 (2*D)-1:-1:D+1]; % Fnodes(1) is a dummy index
  Onode = 2*D;
  ss = 2*D;
end

intra = zeros(ss);
intra(Ops, Onode) = 1;
for d=1:D-1
  if all_Q_to_Qs
    intra(Qnodes(d), Qnodes(d+1:end)) = 1;
  else
    intra(Qnodes(d), Qnodes(d+1)) = 1;
  end
end
for d=D:-1:3
  intra(Fnodes(d), Fnodes(d-1)) = 1;
end
if F1
  intra(Fnodes(2), Fnodes(1)) = 1;
end
if all_Q_to_Qs
  if F1
    intra(Qnodes(1), Fnodes(1:end)) = 1;
  else
    intra(Qnodes(1), Fnodes(2:end)) = 1;
  end
  for d=2:D
    intra(Qnodes(d), Fnodes(d:end)) = 1;
  end
else
  if F1
    intra(Qnodes(1), Fnodes([1 2])) = 1;
  else
    intra(Qnodes(1), Fnodes(2)) = 1;
  end
  for d=2:D-1
    intra(Qnodes(d), Fnodes([d d+1])) = 1;
  end
  intra(Qnodes(D), Fnodes(D)) = 1;
end


inter = zeros(ss);
for d=1:D
  inter(Qnodes(d), Qnodes(d)) = 1;
end
if F1
  inter(Fnodes(1), Qnodes(1)) = 1;
end
for d=2:D
  inter(Fnodes(d), Qnodes([d-1 d])) = 1;
end

if ~F1
  Fnodes = Fnodes(2:end); % strip off dummy -1 term
end
