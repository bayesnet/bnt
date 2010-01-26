
%%%%%%%%%%

function [theta, eta] = mixture_of_experts(q, data, num_iter, theta, eta)
% MIXTURE_OF_EXPERTS Fit a piecewise linear regression model using stochastic gradient descent.
% [theta, eta] = mixture_of_experts(q, data, num_iter)
%
% Inputs:
% q = number of pieces (experts)
% data(l,:) = input example l 
% 
% Outputs:
% theta(i,:) = regression vector for expert i
% eta(i,:) = softmax (gating) params for expert i

[num_cases dim] = size(data);
data = [ones(num_cases,1) data]; % prepend with offset
mu = 0.5; % step size
sigma = 1; % variance of noise

if nargin < 4
  theta = 0.1*rand(q, dim);
  eta = 0.1*rand(q, dim);
end

for t=1:num_iter
  for iter=1:num_cases
    x = data(iter, 1:dim);
    ystar = data(iter, dim+1); % target
    % yhat(i) = E[y | Q=i, x] = prediction of i'th expert
    yhat = theta * x'; 
    % gate_prior(i,:) = Pr(Q=i | x)
    gate_prior = exp(eta * x');
    gate_prior = gate_prior / sum(gate_prior);
    % lik(i) = Pr(y | Q=i, x)
    lik = (1/(sqrt(2*pi)*sigma)) * exp(-(0.5/sigma^2) * ((ystar - yhat) .* (ystar - yhat)));
    % gate_posterior(i,:) = Pr(Q=i | x, y)
    gate_posterior = gate_prior .* lik;
    gate_posterior = gate_posterior / sum(gate_posterior);
    % Update
    eta = eta + mu*(gate_posterior - gate_prior)*x;
    theta = theta + mu*(gate_posterior .* (ystar - yhat))*x;
  end

  if mod(t,100)==0
    fprintf(1, 'iter %d\n', t);
  end

end
fprintf(1, '\n');

