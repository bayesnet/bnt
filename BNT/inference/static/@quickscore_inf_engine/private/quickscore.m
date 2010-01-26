function prob = quickscore(fpos, fneg, inhibit, prior, leak)
% QUICKSCORE Heckerman's algorithm for BN2O networks.
% prob = quickscore(fpos, fneg, inhibit, prior, leak)
% 
% Consider a BN2O (Binary Node 2-layer Noisy-or) network such as QMR with
% dieases on the top and findings on the bottom. (We assume all findings are observed,
% since hidden leaves can be marginalized away.)
% This algorithm takes O(2^|fpos|) time to compute the marginal on all the diseases.
%
% Inputs:
% fpos = the positive findings (a vector of numbers in {1, ..., Nfindings})
% fneg = the negative findings (a vector of numbers in {1, ..., Nfindings})
% inhibit(i,j) = inhibition prob. for finding i, disease j, or 1.0 if j is not a parent.
% prior(j) = prior prob. disease j is ON. We assume prior(off) = 1-prior(on).
% leak(i) = inhibition prob. for the leak node for finding i
%
% Output:
% prob(d) = Pr(disease d = on | ev)
%
% For details, see
% - Heckerman, "A tractable inference algorithm for diagnosing multiple diseases", UAI89.
% - Rish and Dechter, "On the impact of causal independence", UCI tech report, 1998.
%
% Note that this algorithm is numerically unstable, since it adds a large number of positive and
% negative terms and hopes that some of them exactly cancel.
%
% For matlab experts, use 'mex' to compile C_quickscore, which has identical behavior to this function.

[nfindings ndiseases] = size(inhibit);

% make the first disease be always on, for the leak term
Pon = [1 prior(:)'];
Poff = 1-Pon;
Uon = [leak(:) inhibit]; % U(f,d) = Pr(f=0|d=1)
Uoff = [leak(:) ones(nfindings, ndiseases)]; % Uoff(f,d) = Pr(f=0|d=0)
ndiseases = ndiseases + 1;

npos = length(fpos);
post = zeros(ndiseases, 2);
% post(d,1) = alpha Pr(d=off), post(d,2) = alpha Pr(d=m)

FP = length(fpos);
%allbits = logical(dec2bitv(0:(2^FP - 1), FP));
allbits = logical(ind2subv(2*ones(1,FP), 1:(2^FP))-1);

for si=1:2^FP
  bits = allbits(si,:);
  fprime = fpos(bits);
  fmask = zeros(1, nfindings);
  fmask(fneg)=1;
  fmask(fprime)=1;
  fmask = logical(fmask);
  p = 1;
  pterm = zeros(1, ndiseases);
  ptermOff = zeros(1, ndiseases);
  ptermOn = zeros(1, ndiseases);
  for d=1:ndiseases
    ptermOff(d) = prod(Uoff(fmask,d));
    ptermOn(d) = prod(Uon(fmask,d));
    pterm(d) = Poff(d)*ptermOff(d) + Pon(d)*ptermOn(d);
  end
  p = prod(pterm);
  sign = (-1)^(length(fprime));
  for d=1:ndiseases
    myp = p / pterm(d);
    post(d,1) = post(d,1) + sign*(myp * ptermOff(d));
    post(d,2) = post(d,2) + sign*(myp * ptermOn(d));
  end
end

post(:,1) = post(:,1) .* Poff(:);
post(:,2) = post(:,2) .* Pon(:);
post = mk_stochastic(post);
prob = post(2:end,2)'; % skip the leak term


