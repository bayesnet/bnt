% SGS p118
% Try learning the structure using an oracle for the cond indep tests

n = 5;

A = 1; B = 2; C = 3; D = 4; E = 5;

G = zeros(n);
G(A,B)=1;
G(B,[C D]) = 1;
G(C,E)=1;
G(D,E)=1;

k = 2;

pdag = learn_struct_pdag_pc('dsep', n, k, G)




if 0
N = 4; 
dag = zeros(N,N);
C = 1; S = 2; R = 3; W = 4;
dag(C,[R S]) = 1;
dag(R,W) = 1;
dag(S,W)=1;

pdag = learn_struct_pdag_pc('dsep', N, 2, dag)
end
