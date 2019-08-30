function [mu, con] = mixvonMises_Mstep(w, Y, YY, YTY, varargin)
% mixvonMises_Mstep Compute MLEs for mixture of Von Mises given expected sufficient statistics
% function [mu, con] = mixvonMises_Mstep(w, Y, YY, YTY, varargin)
%
% The update for the con parameter only handles the case when ~tied 
% covariance and cov_type ~='s'.(lines 63-74)
% The update for mu only handles ~isempty(clamped_mean)==false and can be
% seen in lines(37-47)
%

[cov_type, tied_cov,  clamped_cov, clamped_mean, cov_prior, other] = ...
    process_options(varargin,...
		    'cov_type', 'full', 'tied_cov', 0,  'clamped_cov', [], 'clamped_mean', [], ...
		    'cov_prior', []);

[Ysz Q] = size(Y);
N = sum(w);
if isempty(cov_prior)
  %cov_prior = zeros(Ysz, Ysz, Q);
  %for q=1:Q
  %  cov_prior(:,:,q) = 0.01*cov(Y(:,q)');
  %end
  cov_prior = repmat(0.01*eye(Ysz,Ysz), [1 1 Q]);
end
%YY = reshape(YY, [Ysz Ysz Q]) + cov_prior; % regularize the scatter matrix
YY = reshape(YY, [Ysz Ysz Q]);

% Set any zero weights to one before dividing
% This is valid because w(i)=0 => Y(:,i)=0, etc
w = w + (w==0);
% This is valid because Y(:,i)=0 => Y(:,:,i)=0, etc
Y = Y + (Y==0);
		    
if ~isempty(clamped_mean)
  mu = clamped_mean;
else
  mu = zeros(Ysz, Q);
  z = zeros(Ysz, Q);
  for i=1:Q
     %take the derivative of sum_t sum_i gamma_t_i(ln w_t+ln
     %2piI_0(k_t)+k_tcos(x_i-mu_0)) with respect to mu which results in the
     %following: arctan(sum_i gamma_t_i sin(x_i)/sum_i gamma_t_i cos(x_i)).
     %Since we already calculated the weighted data previously we can
     %subsitute it in. 
     mu(:,i) = atan2(YY(:,:,i),Y(:,i)); 
  end
end

if ~isempty(clamped_cov)
  con = clamped_cov;
  return;
end

if ~tied_cov
  con = zeros(Ysz,Ysz,Q);
  for i=1:Q
    if cov_type(1) == 's'
      % eqn 17
      error('not complete')
      s2 = (1/Ysz)*( (YTY(i)/w(i)) - mu(:,i)'*mu(:,i) );
      con(:,:,i) = s2 * eye(Ysz);
    else
      % take the derivative of sum_t sum_i gamma_t_i(ln w_t+ln
      %2piI_0(k_t)+k_tcos(x_i-mu_0)) with respect to k. Use the ESS to fill
      %in the blanks and then, find the zeros of the function.
      A = (Y(:,i)*cos(mu(:,i))+YY(:,:,i)*sin(mu(:,i)))/w(i);
      %check that A is within 0 and 1. 
      if (A>=0 && A<1-0.001)
        SS = fzero(@(num) (besseli(1,num)/besseli(0,num))-A,[0,100]);
        %approximation Bannerjee et al.
        %SS = A*(2-A^2)/(1-A^2);
      else
        SS = 1;
      end
      if cov_type(1)=='d'
        SS = diag(diag(SS));
      end
      con(:,:,i) = SS;
    end
  end
else % tied cov
  if cov_type(1) == 's'
    % eqn 19
    s2 = (1/(N*Ysz))*(sum(YTY,2) + sum(diag(mu'*mu) .* w));
    con = s2*eye(Ysz);
  else
    SS = zeros(Ysz, Ysz);
    % eqn 15
    for i=1:Q % probably could vectorize this...
      SS = SS + YY(:,:,i)/N - mu(:,i)*mu(:,i)';
    end
    if cov_type(1) == 'd'
      con = diag(diag(SS));
    else
      con = SS;
    end
  end
end

if tied_cov
  con =  repmat(con, [1 1 Q]);
end
con = con + cov_prior;
