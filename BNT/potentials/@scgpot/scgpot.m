function pot = scgpot(ddom, cheaddom, ctaildom, node_sizes, scgpotc)
% SCGPOT Make a stable CG potential.
% pot = scgpot(ddom, cheaddom, ctaildom, node_sizes, scgpotc)
%
% ddom is discrete nodes contains in the potential
% cheaddom is head nodes constains in the potential
% ctaildom is tail nodes contains in the potential
% node_sizes(i) is the size of the i'th node.
% scgpotc is list of scgcpot objects.

pot.ddom = ddom;
pot.cheaddom = cheaddom;
pot.ctaildom = ctaildom;
pot.domain = myunion(ddom, myunion(cheaddom, ctaildom));
pot.dsizes = node_sizes(pot.ddom);
pot.dsize = prod(node_sizes(pot.ddom));
pot.cheadsizes = node_sizes(pot.cheaddom);
pot.cheadsize = sum(node_sizes(pot.cheaddom));
pot.ctailsizes = node_sizes(pot.ctaildom);
pot.ctailsize = sum(node_sizes(pot.ctaildom));

if nargin < 5
    scgpotc = cell(1, pot.dsize);
    for i=1:pot.dsize
        scgpotc{i} = scgcpot(pot.cheadsize, pot.ctailsize);
    end
end
pot.scgpotc = scgpotc;              

pot = class(pot, 'scgpot');
