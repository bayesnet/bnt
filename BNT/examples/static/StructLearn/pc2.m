% SGS p141 (female orgasm data set)

C = eye(7,7);
C(2,1:1) = [-0.132];
C(3,1:2) = [0.009 -0.136];
C(4,1:3) = [0.22 -0.166 0.403];
C(5,1:4) = [-0.008 0.008 0.598 0.282];
C(6,1:5) = [0.119 -0.076 0.264 0.514 0.176];
C(7,1:6) = [0.118 -0.137 0.368 0.414 0.336 0.338];

n = 7;
for i=1:n
  for j=i+1:n
    C(i,j)=C(j,i);
  end
end

max_fan_in = 4;
nsamples = 281;
alpha = 0.05;
pdag = learn_struct_pdag_pc('cond_indep_fisher_z', n, max_fan_in, C, nsamples, alpha)
