function CPD = set_fields(CPD, varargin)
% SET_PARAMS Set the parameters (fields) for a tabular_CPD object
% CPD = set_params(CPD, name/value pairs)
%
% The following optional arguments can be specified in the form of name/value pairs:
%
% CPT     - the CPT
% prior   - the prior
% clamped - 1 means don't adjust during EM
%
% e.g., CPD = set_params(CPD, 'CPT', 'rnd')

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'CPT', 
    if ischar(args{i+1})
      switch args{i+1}
       case 'unif', CPD.CPT = mk_stochastic(myones(CPD.sizes));
       case 'rnd',  CPD.CPT = mk_stochastic(myrand(CPD.sizes));
       otherwise,   error(['invalid type ' args{i+1}]);       
      end
    elseif isscalarBNT(args{i+1})
      p = args{i+1};
      k = CPD.sizes(end);
      % Bug fix by Hervé BOUTROUILLE 10/1/01
      CPD.CPT = myreshape(sample_dirichlet(p*ones(1,k), prod(CPD.sizes(1:end-1)), CPD.sizes));   
      %CPD.CPT = myreshape(sample_dirichlet(p*ones(1,k), prod(CPD.sizes(1:end-1))), CPD.sizes);
    else
      CPD.CPT = myreshape(args{i+1}, CPD.sizes);
    end
   
   case 'prior',       
    if ischar(args{i+1}) & strcmp(args{i+1}, 'unif')
      CPD.prior = myones(CPD.sizes);
    elseif isscalarBNT(args{i+1})
      CPD.prior = args{i+1} * normalise(myones(CPD.sizes));
    else
      CPD.prior = myreshape(args{i+1}, CPD.sizes);
    end
    
   %case 'clamped',      CPD.clamped = strcmp(args{i+1}, 'yes');
   %case 'clamped',      CPD = set_clamped(CPD, strcmp(args{i+1}, 'yes'));
   case 'clamped',      CPD = set_clamped(CPD, args{i+1});
   
   otherwise,  
    %error(['invalid argument name ' args{i}]);       
  end
end


