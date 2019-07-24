% Let the user draw a square with the mouse,
% and then click on the corners to do a manual segmentation

ss = 6;
Q1 = 1; Q2 = 2; Q3 = 3; obsvel = 6;
CLOCKWISE = 1; ANTICLOCK = 2;
LR = 1; UD = 2; RL = 3; DU = 4;

% repeat this block manually incrementing the sequence number
% and setting ori.
% (since I don't know how to call getmouse as a call-return function).
seq = 4;
%ori = CLOCKWISE
ori = ANTICLOCK;
clear xpos ypos
getmouse
% end block

% manual segmentation with the mouse
startseg(1) = 1;
for i=2:4
  fprintf('click on start of segment %d\n', i);
  [x,y] = ginput(1);
  plot(x,y,'ro')
  d = dist2([xpos; ypos]', [x y]);
  startseg(i) = argmin(d);
end

% plot corners in green 
%ti = first point in (i+1)st segment
t1 = startseg(1); t2 = startseg(2); t3 = startseg(3); t4 = startseg(4); 
plot(xpos(t2), ypos(t2), 'g*')
plot(xpos(t3), ypos(t3), 'g*')
plot(xpos(t4), ypos(t4), 'g*')


xvel = xpos(2:end) - xpos(1:end-1);
yvel = ypos(2:end) - ypos(1:end-1);
speed = [xvel(:)'; yvel(:)'];
pos_data{seq} = [xpos(:)'; ypos(:)'];
vel_data{seq} = [xvel(:)'; yvel(:)'];
T = length(xvel);
Q1label{seq} = num2cell(repmat(ori, 1, T));
Q2label{seq} = zeros(1, T);
if ori == CLOCKWISE
  Q2label{seq}(t1:t2) = LR;
  Q2label{seq}(t2+1:t3) = UD;
  Q2label{seq}(t3+1:t4) = RL;
  Q2label{seq}(t4+1:T) = DU;
else
  Q2label{seq}(t1:t2) = RL;
  Q2label{seq}(t2+1:t3) = UD;
  Q2label{seq}(t3+1:t4) = LR;
  Q2label{seq}(t4+1:T) = DU;
end

% pos_data{seq}(:,t), vel_data{seq}(:,t) Q1label{seq}(t) Q2label{seq}(t)
save 'square4' pos_data vel_data Q1label Q2label

nseq = 4;
cases = cell(1,nseq);
for seq=1:nseq
  T = size(vel_data{seq},2);
  ev = cell(ss,T);
  ev(obsvel,:) = num2cell(vel_data{seq},1);
  ev(Q1,:) = Q1label{seq};
  ev(Q2,:) = num2cell(Q2label{seq});
  cases{seq} = ev;
end
save 'square4_cases' cases 
