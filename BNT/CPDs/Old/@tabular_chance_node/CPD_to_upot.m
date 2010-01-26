function pot = CPD_to_upot(CPD, domain)
% CPD_TO_UPOT Convert a CPD to a utility potential
% pot = CPD_to_upot(CPD, domain)

sz = CPD.size; % mysize(CPD.CPT);
pot = upot(domain, sz, CPD.CPT, 0*myones(sz));
