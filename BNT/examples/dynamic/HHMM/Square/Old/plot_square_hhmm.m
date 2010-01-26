function plot_square_hhmm(ev)
% Plot the square shape implicit in the evidence.
% ev{i,t} is the value of node i in slice t.
% The observed node contains a velocity (delta increment), which is converted
% into a position.
% The Q2 node specifies which model is used; each segment is color-coded
% in the order red, green, blue, black.

Q1 = 1; Q2 = 2; Q3 = 3; F3 = 4; F2 = 5; Onode = 6;

delta = cell2num(ev(Onode,:)); % delta(:,t)
Q2label = cell2num(ev(Q2,:)); 

T = size(delta, 2);
pos = zeros(2,T+1);
clf
hold on
cols = {'r', 'g', 'b', 'k'};
boundary = 0;
coli = 1;
for t=2:T+1
  pos(:,t) = pos(:,t-1) + delta(:,t-1);
  plot(pos(1,t), pos(2,t), sprintf('%c.', cols{coli}));
  if t < T
    boundary = (Q2label(t) ~= Q2label(t-1));
  end
  if boundary
    coli = coli + 1;
    coli = mod(coli-1, length(cols)) + 1;
  end
end

