function bnet = mk_uffe_dbn()

% Make the Uffe DBN from fig 3.4 p55 of my thesis

ss = 4;
intra = zeros(ss,ss);
intra(1,[2 3])=1;
intra(2,3)=1;
intra(3,4)=1;
inter = zeros(ss,ss);
inter(1,1)=1;
inter(4,4)=1;
ns = 2*ones(1,ss);
bnet = mk_dbn(intra, inter, ns);
