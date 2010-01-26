function CPD = update_CPT(CPD)
% Compute the big CPT for an HHMM F node given internal termprob
% function CPD = update_CPT(CPD)

Qsz = CPD.Qsizes(CPD.Q);
Qpsz = prod(CPD.Qsizes(CPD.Qps));

% P(Q(1:d-1), Q(d), F(d+1), F(d))
CPT = zeros(Qpsz, Qsz, 2, 2);
CPT(:,:,1,1) = 1; % if F(d+1)=1, then F(d)=1
CPT(:,:,2,:) = CPD.termprob;

CPD = set_fields(CPD, 'CPT', CPT);          
