function CPD = maximize_params(CPD, temp)
% MAXIMIZE_PARAMS Find ML params of an MLP using Scaled Conjugated Gradient (SCG)
% CPD = maximize_params(CPD, temperature)
% temperature parameter is ignored

if ~adjustable_CPD(CPD), return; end
options = foptions;

% options(1) >= 0 means print an annoying message when the max. num. iter. is reached
if CPD.verbose
  options(1) = 1;
else
  options(1) = -1;
end
%options(1) = CPD.verbose;

options(2) = CPD.wthresh;
options(3) = CPD.llthresh;
options(14) = CPD.max_iter;

dpsz=length(CPD.mlp);

for i=1:dpsz
    mask=[];
    mask=find(CPD.eso_weights(:,:,i)>0);    % for adapting the parameters we use only positive weighted example
    if  ~isempty(mask),
        CPD.mlp{i} = netopt_weighted(CPD.mlp{i}, options, CPD.parent_vals(mask',:), CPD.self_vals(mask',:,i), CPD.eso_weights(mask',:,i), 'scg');
        
        CPD.W1(:,:,i)=CPD.mlp{i}.w1;        % update the parameters matrix
        CPD.b1(i,:)=CPD.mlp{i}.b1;          %
        CPD.W2(:,:,i)=CPD.mlp{i}.w2;        % update the parameters matrix
        CPD.b2(i,:)=CPD.mlp{i}.b2;          %
    end
end
