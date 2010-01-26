% Plot how precision matrix changes over time for KF solution

seed = 0;
rand('state', seed);
randn('state', seed);

[A,B,C,Q,R,Qbig,Rbig,init_x,init_V,robot_block,landmark_block,...
	  true_landmark_pos, true_robot_pos, true_data_assoc, ...
	  obs_rel_pos, ctrl_signal] = mk_linear_slam(...
	      'nlandmarks', 6, 'T', 12, 'ctrl', 'leftright', 'data-assoc', 'cycle');

figure(1); clf
hold on
for i=1:nlandmarks
  %text(true_landmark_pos(1,i), true_landmark_pos(2,i), sprintf('L%d',i));
  plot(true_landmark_pos(1,i), true_landmark_pos(2,i), '*')
end
hold off


[x, V] = kalman_filter(obs_rel_pos, A, C, Qbig, Rbig, init_x, init_V, ...
				     'model', true_data_assoc, 'u', ctrl_signal, 'B', B);

est_robot_pos = x(robot_block, :);
est_robot_pos_cov = V(robot_block, robot_block, :);

for i=1:nlandmarks
  bi = landmark_block(:,i);
  est_landmark_pos(:,i) = x(bi, T);
  est_landmark_pos_cov(:,:,i) = V(bi, bi, T);
end


if 0
figure(1); hold on
for i=1:nlandmarks
  h=plotgauss2d(est_landmark_pos(:,i), est_landmark_pos_cov(:,:,i));
  set(h, 'color', 'r')
end
hold off

hold on
for t=1:T
  h=plotgauss2d(est_robot_pos(:,t), est_robot_pos_cov(:,:,t));
  set(h,'color','r')
  h=text(est_robot_pos(1,t), est_robot_pos(2,2), sprintf('R%d', t));
  set(h,'color','r')
end
hold off
end


P = zeros(size(V));
for t=1:T
  P(:,:,t) = inv(V(:,:,t));
end

if 0
  figure(2)
  for t=1:T
    subplot(T/2,2,t)
    imagesc(P(1:2:end,1:2:end, t))
    colorbar
  end
else
  figure(2)
  for t=1:T
    subplot(T/2,2,t)
    imagesc(V(1:2:end,1:2:end, t))
    colorbar
  end
end

% marginalize out robot position and then check structure
bi = landmark_block(:);
V = V(bi,bi,T); 
P = inv(V);
P(1:2:end,1:2:end)
