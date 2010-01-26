function plot_square_hhmm(ev)
% Plot the square shape implicit in the evidence.
% ev{i,t} is the value of node i in slice t.
% The observed node contains a velocity (delta increment), which is converted
% into a position.
% The Q2 node specifies which model is used, and hence which color
% to use: 1=red, 2=green, 3=blue, 4=black.

Q1 = 1; Q2 = 2; Q3 = 3; F3 = 4; F2 = 5; Onode = 6;

delta = cell2num(ev(Onode,:)); % delta(:,t)
Q2label = cell2num(ev(Q2,:)); 

T = size(delta, 2);
pos = zeros(2,T+1);
hold on
cols = {'r', 'g', 'b', 'k'};
for t=2:T+1
  pos(:,t) = pos(:,t-1) + delta(:,t-1);
  plot(pos(1,t), pos(2,t), sprintf('%c.', cols{Q2label(t-1)}));
  if (t==2)
    text(pos(1,t-1),pos(2,t-1),sprintf('%d',t))
  elseif (mod(t,20)==0)
    text(pos(1,t),pos(2,t),sprintf('%d',t))
  end
end

