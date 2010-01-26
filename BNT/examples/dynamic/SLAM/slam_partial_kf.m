% See how well partial Kalman filter updates work

seed = 0;
rand('state', seed);
randn('state', seed);
nlandmarks = 6;
T = 12;

[A,B,C,Q,R,Qbig,Rbig,init_x,init_V,robot_block,landmark_block,...
	  true_landmark_pos, true_robot_pos, true_data_assoc, ...
	  obs_rel_pos, ctrl_signal] = mk_linear_slam(...
	      'nlandmarks', nlandmarks, 'T', T, 'ctrl', 'leftright', 'data-assoc', 'cycle');

% exact
[xe, Ve] = kalman_filter(obs_rel_pos, A, C, Qbig, Rbig, init_x, init_V, ...
				     'model', true_data_assoc, 'u', ctrl_signal, 'B', B);


% approx
%k = nlandmarks-1; % exact
k = 3;
ndx = {};
for t=1:T
  landmarks = unique(true_data_assoc(t:-1:max(t-k,1)));
  tmp = [landmark_block(:, landmarks) robot_block'];
  ndx{t} = tmp(:);
end

[xa, Va] = kalman_filter(obs_rel_pos, A, C, Qbig, Rbig, init_x, init_V, ...
			 'model', true_data_assoc, 'u', ctrl_signal, 'B', B, ...
		       'ndx', ndx);



nrows = 10;
stepsize = T/(2*nrows);
ts = 1:stepsize:T;

if 1 % plot
  
clim = [0 max(max(Va(:,:,end)))];

figure(2)
if 0
  imagesc(Ve(1:2:end,1:2:end, T))
  clim = get(gca,'clim');
else
  i = 1;
  for t=ts(:)'
    subplot(nrows,2,i)
    i = i + 1;
    imagesc(Ve(1:2:end,1:2:end, t))
    set(gca, 'clim', clim)
    colorbar
  end
end
suptitle('exact')


figure(3)
if 0
  imagesc(Va(1:2:end,1:2:end, T))
  set(gca,'clim', clim)
else
  i = 1;
  for t=ts(:)'
    subplot(nrows,2,i)
    i = i+1;
    imagesc(Va(1:2:end,1:2:end, t))
    set(gca, 'clim', clim)
    colorbar
  end
end
suptitle('approx')


figure(4)
i = 1;
for t=ts(:)'
  subplot(nrows,2,i)
  i = i+1;
  Vd = Va(1:2:end,1:2:end, t) - Ve(1:2:end,1:2:end,t);
  imagesc(Vd)
  set(gca, 'clim', clim)
  colorbar
end
suptitle('diff')

end % all plot


for t=1:T
  %err(t)=rms(xa(:,t), xe(:,t));
  err(t)=rms(xa(1:end-2,t), xe(1:end-2,t)); % exclude robot
end
figure(5);plot(err)
title('rms mean pos')


for t=1:T
  i = 1:2*nlandmarks;
  denom = Ve(i,i,t) + (Ve(i,i,t)==0);
  Vd =(Va(i,i,t)-Ve(i,i,t)) ./ denom;
  Verr(t) = max(Vd(:));
end
figure(6); plot(Verr)
title('max relative Verr')
