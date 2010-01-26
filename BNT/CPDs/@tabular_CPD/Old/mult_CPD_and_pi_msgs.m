function T = mult_CPD_and_pi_msgs(CPD, n, ps, msgs, except)
% MULT_CPD_AND_PI_MSGS Multiply the CPD and all the pi messages from parents, perhaps excepting one
% T = mult_CPD_and_pi_msgs(CPD, n, ps, msgs, except)

if nargin < 5, except = -1; end

dom = [ps n];
%ns = sparse(1, max(dom));
ns = zeros(1, max(dom));
ns(dom) = mysize(CPD.CPT);
T = dpot(dom, ns(dom), CPD.CPT);
for i=1:length(ps)
  p = ps(i);
  if p ~= except
    T = multiply_by_pot(T, dpot(p, ns(p), msgs{n}.pi_from_parent{i}.T)); 
  end
end         
