% Compare Kalman smoother with loopy

seed = 0;
rand('state', seed);
randn('state', seed);
nlandmarks = 6;
T = 12;

[A,B,C,Q,R,Qbig,Rbig,init_x,init_V,robot_block,landmark_block,...
	  true_landmark_pos, true_robot_pos, true_data_assoc, ...
	  obs_rel_pos, ctrl_signal] = mk_linear_slam(...
	      'nlandmarks', nlandmarks, 'T', T, 'ctrl', 'leftright', 'data-assoc', 'cycle');

[xsmooth, Vsmooth] = kalman_smoother(obs_rel_pos, A, C, Qbig, Rbig, init_x, init_V, ...
				     'model', true_data_assoc, 'u', ctrl_signal, 'B', B);

est_robot_pos = xsmooth(robot_block, :);
est_robot_pos_cov = Vsmooth(robot_block, robot_block, :);

for i=1:nlandmarks
  bi = landmark_block(:,i);
  est_landmark_pos(:,i) = xsmooth(bi, T);
  est_landmark_pos_cov(:,:,i) = Vsmooth(bi, bi, T);
end


if 1
  [bnet, Unode, Snode, Lnodes, Rnode, Ynode, Lsnode] = ...
      mk_gmux_robot_dbn(nlandmarks, Q, R, init_x, init_V, robot_block, landmark_block);
  engine = pearl_unrolled_dbn_inf_engine(bnet, 'max_iter', 50, 'filename', ...
					 '/home/eecs/murphyk/matlab/loopyslam.txt');
else
  [bnet, Unode, Snode, Lnodes, Rnode, Ynode] = ...
      mk_gmux2_robot_dbn(nlandmarks, Q, R, init_x, init_V, robot_block, landmark_block);
  engine = jtree_dbn_inf_engine(bnet);
end

nnodes = bnet.nnodes_per_slice;
evidence = cell(nnodes, T);
evidence(Ynode, :) = num2cell(obs_rel_pos, 1);
evidence(Unode, :) = num2cell(ctrl_signal, 1);
evidence(Snode, :) = num2cell(true_data_assoc);

[engine, ll, niter] = enter_evidence(engine, evidence);
niter

loopy_est_robot_pos = zeros(2, T);
for t=1:T
  m = marginal_nodes(engine, Rnode, t);
  loopy_est_robot_pos(:,t) = m.mu;
end

for i=1:nlandmarks
  m = marginal_nodes(engine, Lnodes(i), T);
  loopy_est_landmark_pos(:,i) = m.mu;
  loopy_est_landmark_pos_cov(:,:,i) = m.Sigma;
end


