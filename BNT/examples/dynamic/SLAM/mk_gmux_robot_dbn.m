function [bnet, Unode, Snode, Lnodes, Rnode, Ynode, Lsnode] = ...
    mk_gmux_robot_dbn(nlandmarks, Q, R, init_x, init_V, robot_block, landmark_block)

% Make DBN

% S
% | L1 -------> L1'
% |  | L2 ----------> L2'
% \  | /
%  v v v
%    Ls
%    |
%    v
%    Y
%    ^
%    |
%    R ------->  R'
%    ^      
%    |      
%    U      
%
%
% S is a switch, Ls is a deterministic gmux, Y = Ls-R,
% R(t+1) = R(t) + U(t+1), L(t+1) = L(t)


% number nodes topologically
Snode = 1;
Lnodes = 2:nlandmarks+1;
Lsnode = nlandmarks+2;
Unode = nlandmarks+3;
Rnode = nlandmarks+4;
Ynode = nlandmarks+5;

nnodes = nlandmarks+5; 
intra = zeros(nnodes, nnodes);
intra([Snode Lnodes], Lsnode) =1;
intra(Unode,Rnode)=1;
intra([Rnode Lsnode], Ynode)=1;

inter = zeros(nnodes, nnodes);
inter(Rnode, Rnode)=1;
for i=1:nlandmarks
  inter(Lnodes(i), Lnodes(i))=1;
end

Lsz = 2; % (x y) posn of landmark
Rsz = 2; % (x y) posn of robot
Ysz = 2; % relative distance
Usz = 2; % (dx dy) ctrl
Ssz = nlandmarks; % can switch between any landmark

ns = zeros(1,nnodes);
ns(Snode) = Ssz;
ns(Lnodes) = Lsz;
ns(Lsnode) = Lsz;
ns(Ynode) = Ysz;
ns(Rnode) = Rsz;
ns(Ynode) = Usz;
ns(Unode) = Usz;

bnet = mk_dbn(intra, inter, ns, 'discrete', Snode, 'observed', [Snode Ynode Unode]);


bnet.CPD{Snode} = root_CPD(bnet, Snode); % always observed
bnet.CPD{Unode} = root_CPD(bnet, Unode); % always observed
for i=1:nlandmarks
  bi = landmark_block(:,i);
  bnet.CPD{Lnodes(i)} = gaussian_CPD(bnet, Lnodes(i), 'mean', init_x(bi), 'cov', init_V(bi,bi));
end
bi = robot_block;
bnet.CPD{Rnode} = gaussian_CPD(bnet, Rnode, 'mean', init_x(bi), 'cov', init_V(bi,bi), 'weights', eye(2));
bnet.CPD{Lsnode} = gmux_CPD(bnet, Lsnode, 'cov', repmat(zeros(Lsz,Lsz), [1 1 nlandmarks]), ...
			    'weights', repmat(eye(Lsz,Lsz), [1 1 nlandmarks]));
W = [eye(2) -eye(2)]; % Y = Ls - R, where Ls is the lower-numbered parent
bnet.CPD{Ynode} = gaussian_CPD(bnet, Ynode, 'mean', zeros(Ysz,1), 'cov', R, 'weights', W);

% slice 2
eclass = bnet.equiv_class;
W = [eye(2) eye(2)]; % R(t) = R(t-1) + U(t), where R(t-1) is the lower-numbered parent
bnet.CPD{eclass(Rnode,2)} = gaussian_CPD(bnet, Rnode+nnodes, 'mean', zeros(Rsz,1), 'cov', Q, 'weights', W);
for i=1:nlandmarks
  bnet.CPD{eclass(Lnodes(i), 2)} = gaussian_CPD(bnet, Lnodes(i)+nnodes, 'mean', zeros(2,1), ...
						   'cov', zeros(2,2), 'weights', eye(2));
end
