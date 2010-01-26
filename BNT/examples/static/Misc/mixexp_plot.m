function plot_mixexp(theta, eta, data)
% PLOT_MIXEXP  Plot the results for a piecewise linear regression model
% plot_mixexp(theta, eta, data)
% 
% data(l,:) = [x y] for example l
% theta(i,:) = regression vector for expert i
% eta(i,:) = softmax (gating) params for expert i

numexp = size(theta, 1);

mn = min(data);
mx = max(data);
xa = mn(1):0.01:mx(1);
x = [ones(length(xa),1) xa'];
% pr(i,l) = posterior probability of expert i on example l
pr = exp(eta * x');
pr = pr ./ (ones(numexp,1) * sum(pr));
% y(i,l) = prediction of expert i for example l
y = theta * x';
% yg(l) = weighted prediction  for example l
yg = sum(y .* pr)';

subplot(3,2,1);
plot(xa, y(1,:));
title('expert 1');

subplot(3,2,2);
plot(xa, y(2,:));
title('expert 2');

subplot(3,2,3);
plot(xa, pr(1,:));
title('gating 1');

subplot(3,2,4);
plot(xa, pr(2,:));
title('gating 2');

subplot(3,2,5);
plot(xa, yg);
axis([-1 1 -1 2])
title('prediction');

subplot(3,2,6);
title('data');
hold on
plot(data(:,1), data(:,2), '+');
hold off

