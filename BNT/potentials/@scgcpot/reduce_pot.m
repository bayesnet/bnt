function [reduced_pot,successful] = reduce(pot,tailnodes)
% Executes the reduce operation defined in 
% Stable Local Computation with Conditional Gaussian Distributions
% Steffen L. Lauritzen
% Frank Jensen
% September 1999
% The potential pot is reduced if B contains any zero columns 
% The test are restricted to the positions in tailnodes.
% Any columns successfully deleted are entered in the array successful

if nargin < 2
    tailnodes = 1:pot.ctailsize;
end

successful = [];

% Look for all columns beeing equal to zero
for i = tailnodes
    if ~any(pot.B(:,i))
        successful = [successful i]; 
    end
end

remain = mysetdiff(1:pot.ctailsize,successful);

% Erase the zero-columns and decrease the tailsize
pot.B = pot.B(:,remain);
pot.ctailsize = pot.ctailsize - length(successful);

% Return the reduced potential
reduced_pot = pot;
  
