function fh=hme_reg_plot(net, nodes_info, train_data, test_data)
% 
% Use this function ONLY when the input dimension is 1
% and the problem is a regression one.
% We assume that each row of 'train_data' & 'test_data' is an example.
%
% ----------------------------------------------------------------------------------------------------
% -> pierpaolo_b@hotmail.com   or   -> pampo@interfree.it
% ----------------------------------------------------------------------------------------------------

fh=figure('Name','HME based regression', 'MenuBar', 'none', 'NumberTitle', 'off');

mn_x_train = round(min(train_data(:,1)));
mx_x_train = round(max(train_data(:,1)));     
x_train = mn_x_train(1):0.01:mx_x_train(1);
Z_train=fhme(net, nodes_info, x_train',size(x_train,2));      % forward propagation trougth the HME

if nargin==4,
    subplot(2,1,1);
    mn_x_test = round(min(test_data(:,1)));
    mx_x_test = round(max(test_data(:,1)));
    x_test = mn_x_test(1):0.01:mx_x_test(1);
    Z_test=fhme(net, nodes_info, x_test',size(x_test,2));      % forward propagation trougth the HME
end

hold on;
set(gca, 'Box', 'on');
plot(x_train', Z_train, 'r');
plot(train_data(:,1),train_data(:,2),'+k');
title('Training set and prediction');
hold off

if nargin==4,
    subplot(2,1,2);
    hold on;
    set(gca, 'Box', 'on');
    plot(x_train', Z_train, 'r');
    if size(test_data,2)==2,
        plot(test_data(:,1),test_data(:,2),'+k');
    end
    title('Test set and prediction');
    hold off
end