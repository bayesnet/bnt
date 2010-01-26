function [W, b] = extract_params(CPD)

% W(X,Y,Q), b(Y,Q)  where Y = ns(self), X = ns(cps), Q = prod(ns(dps))

glimsz = prod(CPD.sizes(CPD.dpndx));
ss = CPD.sizes(end);
cpsz       = sum(CPD.sizes(CPD.cpndx));
dp_as_cpsz = sum(CPD.sizes(CPD.dps_as_cps.ndx));
W = zeros(dp_as_cpsz + cpsz, ss, glimsz);
b = zeros(ss, glimsz);

for i=1:glimsz
  W(:,:,i) = CPD.glim{i}.w1;
  b(:,i) = CPD.glim{i}.b1(:);
end

W = myreshape(W, [dp_as_cpsz + cpsz ss CPD.sizes(CPD.dpndx)]);
b = myreshape(b, [ss CPD.sizes(CPD.dpndx)]);
