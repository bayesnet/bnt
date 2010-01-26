function [reduced_pot,successful] = reduce_pot(pot,tailnodes)
% Executes the reduce operation defined in
% Stable Local Computation with Conditional Gaussian Distributions
% Steffen L. Lauritzen
% Frank Jensen
% September 1999
% The potential pot is reduced if B contains any zero columns
% The test are restricted to the positions in tailnodes.
% Any columns successfully deleted are entered in the array successful
if nargin < 2
    tailnodes = pot.ctaildom;
end

successful = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keep track of remaining tailnodes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rem_tailnodes = pot.ctaildom;
for i = tailnodes
    pos = find(i==rem_tailnodes);
    successful_red = [pos];
    red_scgcpot = cell(1,pot.dsize);
    j = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Test whether all components of pot.scgpotc can be reduced %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while ((j <= pot.dsize) & ~isempty(successful_red))
        [cpot,successful_red] = reduce_pot(pot.scgpotc{j},pos);
        red_scgcpot{j} = cpot;
        j = j + 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If i is a reducible tailnode, then reduce the potential %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(successful_red)
        successful = [successful i];
        pot.scgpotc = red_scgcpot;
        rem_tailnodes = mysetdiff(rem_tailnodes,i);
    end;
end

pot.ctaildom = rem_tailnodes;
positions = find_equiv_posns(rem_tailnodes,pot.ctaildom);
pot.ctailsizes = pot.ctailsizes(positions);
pot.ctailsize = sum(pot.ctailsizes);
pot.domain = mysetdiff(pot.domain,successful);
reduced_pot = pot;





