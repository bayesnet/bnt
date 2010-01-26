function pot = cgpot(ddom, cdom, node_sizes, can, mom, subtype)
% CPOT Make a canonical CG potential.
% function pot = cgpot(ddom, cdom, node_sizes, can, mom, subtype)
%
% node_sizes(i) is the size of the i'th node.
% can and mom default to 0s.
% subtype defaults to 'c'.

if nargin < 6, subtype = 'c'; end

pot.ddom = ddom;
pot.cdom = cdom;
node_sizes = node_sizes(:)'; % row vectors print better
pot.domain = myunion(ddom, cdom);
pot.dsizes = node_sizes(pot.ddom);
pot.dsize = prod(node_sizes(pot.ddom));
pot.csizes = node_sizes(pot.cdom);
pot.csize = sum(node_sizes(pot.cdom));
pot.subtype = subtype;

if nargin < 4
  can = cell(1, pot.dsize);
  for i=1:pot.dsize
    can{i} = cpot(cdom, node_sizes(cdom));
  end
end
pot.can = can;              

if nargin < 5
  mom = cell(1, pot.dsize);
  for i=1:pot.dsize
    mom{i} = mpot(cdom, node_sizes(cdom));
  end
end
pot.mom = mom;

pot = class(pot, 'cgpot');

