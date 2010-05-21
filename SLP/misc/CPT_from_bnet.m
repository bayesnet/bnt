function CPT = CPT_from_bnet(bnet,v)
% Export CPT from bnet
%   CPT = CPT_from_bnet(bnet,v)
%
% optional : v~=0 --> verbose mode

if nargin<2, v=0; end
N = size(bnet.dag,1);
CPT = cell(1,N);
for j=1:N
  CPD=struct(bnet.CPD{j});
  %counts{j}=CPD.counts;
  CPT{j}=CPD.CPT;
end
if v, celldisp(CPT);end
%nsamples = CPD.nsamples;
