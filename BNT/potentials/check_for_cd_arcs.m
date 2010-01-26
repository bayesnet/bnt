function check_for_cd_arcs(onodes, cnodes, dag)
% CHECK_FOR_CD_ARCS Raise an error if there are any C->D links where the C node is hidden.
% check_for_cd_arcs(onodes, cnodes, dag)
%
% We cannot convert the logistic/softmax function (C->D CPD) to a Gaussian potential
% unless we use the variational approximation discussed in 
% "A variational approximation for Bayesian networks with discrete and continuous latent
% variables", K. Murphy, UAI 1999.

n = length(dag);
hnodes = mysetdiff(1:n, onodes);
chid = myintersect(cnodes, hnodes);
dnodes = mysetdiff(1:n, cnodes);
for i=chid(:)'
  dcs = myintersect(children(dag, i), dnodes);
  if ~isempty(dcs)
    error(['hidden cts node ' num2str(i) ' has a discrete child']);
  end
end




