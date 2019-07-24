function CPD = reset_ess(CPD)
% RESET_ESS Reset the Expected Sufficient Statistics of a hhmm2 Q node.
% CPD = reset_ess(CPD)

domsz = CPD.dom_sz;
domsz(CPD.Fself_ndx) = 1;
domsz(CPD.Fbelow_ndx) = 1;
Qdom_sz = domsz;
Qdom_sz(Qdom_sz==1)=[]; % get rid of dimensions of size 1

CPD.start_counts = zeros(Qdom_sz);
CPD.trans_counts = zeros(Qdom_sz);
