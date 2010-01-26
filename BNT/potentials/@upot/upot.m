function pot = upot(domain, sizes, p, u)
% UPOT Make a discrete utility potential.
% pot = upot(domain, sizes, p, u)
%
% sizes(i) is the size of the i'th domain element.
% p defaults to all 1s, u defaults to all 0s.

if nargin < 3, p = myones(sizes); end
if nargin < 4, u = 0*myones(sizes); end

pot.domain = domain;
pot.p = myreshape(p, sizes);
pot.u = myreshape(u, sizes);
pot.sizes = sizes(:)';
pot = class(pot, 'upot');
