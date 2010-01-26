function pot = dpot(domain, sizes, T)
% DPOT Make a discrete (sparse) potential.
% pot = dpot(domain, sizes, T, spar)
%
% sizes(i) is the size of the i'th domain element.
% T defaults to all 1s.

%assert(length(sizes) == length(domain));

pot.domain = domain(:)'; % so we can see it when we display
if nargin < 3
  pot.T = myones(sizes);
  %pot.T = ones(1,prod(sizes)); % 1D vector
else 
   if isempty(T)
      pot.T = [];
   else
      if issparse(T)
         pot.T = T;   
      else
         pot.T = myreshape(T, sizes);  
      end
   end
end
pot.sizes = sizes(:)';
pot = class(pot, 'dpot');
