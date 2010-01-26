function fh=hme_class_plot(net, nodes_info, train_data, test_data)
% 
% Use this function ONLY when the input dimension is 2
% and the problem is a classification one.
% We assume that each row of 'train_data' & 'test_data' is an example.
%
%------Line Spec------------------------------------------------------------------------
%
% LineWidth       - specifies the width (in points) of the line
% MarkerEdgeColor - specifies the color of the marker or the edge color
%                   forfilled markers (circle, square, diamond, pentagram, hexagram, and the
%                   four triangles).
% MarkerFaceColor - specifies the color of the face of filled markers.
% MarkerSize      - specifies the size of the marker in points.
%
% Example
% -------
% plot(t,sin(2*t),'-mo',...
%                'LineWidth',2,...
%                'MarkerEdgeColor','k',...                          % 'k'=black
%                'MarkerFaceColor',[.49 1 .63],...                  % RGB color
%                'MarkerSize',12)
%----------------------------------------------------------------------------------------

class_num=nodes_info(2,end);
mn_x = round(min(train_data(:,1)));     mx_x = round(max(train_data(:,1)));
mn_y = round(min(train_data(:,2)));     mx_y = round(max(train_data(:,2)));
if nargin==4,
    mn_x = round(min([train_data(:,1); test_data(:,1)]));     
    mx_x = round(max([train_data(:,1); test_data(:,1)]));
    mn_y = round(min([train_data(:,2); test_data(:,2)]));
    mx_y = round(max([train_data(:,1); test_data(:,2)]));
end
x = mn_x(1)-1:0.2:mx_x(1)+1;
y = mn_y(1)-1:0.2:mx_y(1)+1;
[X, Y] = meshgrid(x,y);
X = X(:); 
Y = Y(:);
num_g=size(X,1);
griglia = [X Y];
rand('state',1);
if class_num<=6,
    colors=['r'; 'g'; 'b'; 'c'; 'm'; 'y'];
else
    colors=rand(class_num, 3);  % each row is an RGB color
end
fh = figure('Name','Data & decision boundaries', 'MenuBar', 'none', 'NumberTitle', 'off');
ms=5;           % Marker Size
if nargin==4,
%    ms=4;       % Marker Size
    subplot(1,2,1);
end
% Plot of train_set -------------------------------------------------------------------------
axis([mn_x-1 mx_x+1 mn_y-1 mx_y+1]);
set(gca, 'Box', 'on');
c_max_train = max(train_data(:,3));
hold on
for m=1:c_max_train,
    app_x=train_data(:,1);
    app_y=train_data(:,2);
    thisX=app_x(train_data(:,3)==m);
    thisY=app_y(train_data(:,3)==m);
    if class_num<=6,
       str_col=[];
       str_col=['o', colors(m,:)];
       plot(thisX, thisY, str_col, 'MarkerSize', ms);
    else
        plot(thisX, thisY, 'o',...
            'LineWidth', 1,...            
            'MarkerEdgeColor', colors(m,:), 'MarkerSize', ms)
    end
end
%---hmefwd_generale(net,data,ndata)-----------------------------------------------------------
Z=fhme(net, nodes_info, griglia, num_g);      % forward propagation trougth the HME
%---------------------------------------------------------------------------------------------
[foo , class] = max(Z'); % 0/1 loss function => we assume that the true class is the one with the
                         % maximum posterior prob.
class = class';
for m = 1:class_num,
  thisX=[]; thisY=[];
  thisX = X(class == m);
  thisY = Y(class == m);
  if class_num<=6,
      str_col=[];
      str_col=['d', colors(m,:)];
      h=plot(thisX, thisY, str_col);      
  else
      h = plot(thisX, thisY, 'd',...
          'MarkerEdgeColor',colors(m,:),...
          'MarkerFaceColor','w');
  end
  set(h, 'MarkerSize', 4);
end
title('Training set and Decision Boundaries (0/1 loss)')
hold off

% Plot of test_set --------------------------------------------------------------------------
if nargin==4,
    subplot(1,2,2);
    axis([mn_x-1 mx_x+1 mn_y-1 mx_y+1]);
    set(gca, 'Box', 'on');
    hold on     
    if size(test_data,2)==3,  % we know the classification of the test set examples
        c_max_test = max(test_data(:,3));
        for m=1:c_max_test,
            app_x=test_data(:,1);
            app_y=test_data(:,2);
            thisX=app_x(test_data(:,3)==m);
            thisY=app_y(test_data(:,3)==m);
            if class_num<=6,
                str_col=[];
                str_col=['o', colors(m,:)];
                plot(thisX, thisY, str_col, 'MarkerSize', ms);
            else
                plot(thisX, thisY, 'o',...
                     'LineWidth', 1,...
                     'MarkerEdgeColor', colors(m,:),...
                     'MarkerSize',ms);
            end
        end
    else
        plot(test_data(:,1), test_data(:,2), 'ko',...
            'MarkerSize', ms);
    end
    for m = 1:class_num,
        thisX=[]; thisY=[];
        thisX = X(class == m);
        thisY = Y(class == m);
        if class_num<=6,
          str_col=[];
          str_col=['d', colors(m,:)];
          h=plot(thisX, thisY, str_col);  
        else
           h = plot(thisX, thisY, 'd',...
                    'MarkerEdgeColor', colors(m,:),...
                    'MarkerFaceColor','w');
        end
        set(h, 'MarkerSize', 4);
    end
    title('Test set and Decision Boundaries (0/1 loss)')
    hold off
end