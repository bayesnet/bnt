function CPD = update_CPT(CPD)
% Compute the big CPT for an HHMM F node given internal termprob
% function CPD = update_CPT(CPD)

Qsz = CPD.Qsz;
Qpsz = CPD.Qpsz;

% CPT(Qpsz, Q, Fbelow, Fself)
CPT = zeros(Qpsz, Qsz, 2, 2);
CPT(:,:,1,1) = 1; % if Fbelow=1 (off), then Fself=1 (off)
CPT(:,:,2,:) = CPD.termprob;

CPD = set_fields(CPD, 'CPT', CPT);          
