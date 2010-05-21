function [Phase_3, Phase_2, Phase_1, UPhase_3] = learn_struct_bnpc(data,node_sizes,epsilon,star)
% G = learn_struct_bnpc(Data,node_sizes,epsilon,star)
%
% Data(i,m) is node i in case m.
% node_sizes and epsilon are optionnals.
% star = 0 to use try_to_separate_B instead of try_to_separate_B_star
%
% see "Learning bayesian Networks from Data: A Efficient Approach Based on Information Theorie"
%     Jie Cheng, David Bell and Weird Liu.
%
% Things to do : rewrite function orient_edges !
%                ! sometimes it causes crashes !
%
% V0.91 : 18 sept 2003 (olivier.francois@insa-rouen.fr)

verbose=1;
%if nargin < 5, mwst=0; end
if nargin < 4, star=1; end
if nargin < 3, epsilon=0.05; end
if nargin < 2, node_sizes=max(data'); end

if verbose
  fprintf('================== phase I : \n');
end
tmp1=cputime;
[Phase_1 II JJ score_mat score_mat2] = phaseI(data, node_sizes, epsilon);
tmp1=cputime-tmp1;

if verbose
  fprintf('Execution time : %2.5f\n',tmp1);
  fprintf('\n================== phase II : \n');
end
tmp1=cputime;
Phase_2 = phaseII(Phase_1, data, node_sizes, epsilon, II, JJ, score_mat);
tmp1=cputime-tmp1;

if verbose
  fprintf('Execution time : %2.5f\n',tmp1);
  fprintf('\n================== phase III : \n');
end
tmp1=cputime;
[Phase_3 UPhase_3]  = phaseIII(Phase_2, data, node_sizes, epsilon, score_mat2, star);
tmp1=cputime-tmp1;
if verbose
  fprintf('Execution time : %2.5f\n',tmp1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G, II, JJ, score_mat, sc2] = phaseI(data,node_sizes,alpha)
% [G, II , JJ, score_mat, s2] = phaseI(data,node_sizes,epsilon)
%
% G is an acyclic graph
% [II JJ] is the list of important edges not processed in phase I (for phase II)
% score_mat is the mutual information score matrix
%
% data(i,m) is node i in case m.
% alpha is the significant level for CI tests ( default=0.05 ).
% node_sizes is the vector of sizes ( default=max(data') ).
%
% see "Learning bayesian Networks from Data: A Efficient Approach Based on Information Theorie"
%     Jie Cheng, David Bell and Weird Liu.

% 0.
if nargin < 3, alpha=0.05; end
if nargin < 2, node_sizes=max(data'); end
[N m] = size(data);
score_mat = zeros(N);
edges=0;

% 1.
G = zeros(N);
L=[];

% 2. Use of Chi2 instead of MI ... allow using a confidence level alpha instead of an arbitrary epsilon
for i=1:(N-1)
  for j=(i+1):N
    [I score_mat(i,j)] = cond_indep_chisquare(i,j,[],data,'LRT',alpha,node_sizes);
  end
end
sc2=score_mat;


[tmp ordre]=sort(-score_mat(:));
ordre2=ordre(find(-tmp>alpha));
[II JJ]=ind2sub([N N],ordre2);

pointer=1 ;
fini=length(II);

% 3.
edges=2;
for pointer=1:min(2,fini),
  %fprintf('%d-%d\n',II(pointer),JJ(pointer));
  G(II(pointer),JJ(pointer))=1;
  G(JJ(pointer),II(pointer))=1;
  score_mat(II(pointer),JJ(pointer))=-inf;
end

pointer=min(2,fini);
arret=0;

while pointer<fini & ~arret
  % 4.
  pointer=pointer+1;
  C = ~reachability_graph(G);
  if C(II(pointer),JJ(pointer))
    %fprintf('%d-%d\n',II(pointer),JJ(pointer));
    G(II(pointer),JJ(pointer))=1;
    G(JJ(pointer),II(pointer))=1;
    score_mat(II(pointer),JJ(pointer))=-inf;
    edges=edges+1;
    if edges==N-1
      arret=1;
    end
  end

  % 5.
end

[tmp ordre]=sort(-score_mat(:));
ordre2=ordre(find(-tmp>alpha));
[II JJ]=ind2sub([N N],ordre2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = phaseII(G1,data,node_sizes,alpha,II, JJ, score_mat)
% G = phaseII(G1,data,node_sizes,epsilon,II,JJ, score_mat)
%
% G1, II, JJ,score_mat are given by phaseI.
% data(i,m) is node i in case m.
% node_sizes is the vector of sizes ( default=max(data') ).
% alpha is the significant level for CI tests ( default=0.05 ).
%
% see "Learning bayesian Networks from Data: A Efficient Approach Based on Information Theorie"
%     Jie Cheng, David Bell and Weird Liu.

st{1}='added';
st{2}='';
% 0.
[N m] = size(data);

G=G1;

% 6.
II=II(end:-1:1);
JJ=JJ(end:-1:1);

pointer=length(II);

while pointer>0
  % 7.
  trysep = try_to_separate_A(G,II(pointer),JJ(pointer),data,alpha,node_sizes);
  %fprintf('Try to separate %d and %d : %s\n',II(pointer),JJ(pointer),st{trysep+1});
  if ~trysep
    G(II(pointer),JJ(pointer))=1;
    G(JJ(pointer),II(pointer))=1;
    %fprintf('%d-%d\n',II(pointer),JJ(pointer));
  end
  % 8.
  pointer=pointer-1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G,U] = phaseIII(G1,data,node_sizes,alpha,s2,star)
% [G] = phaseIII(G1,data,node_sizes,alpha,s2,star)
%
% data(i,m) is node i in case m.
% if star~=0, use try_to_separate_B_star instead of try_to_separate_B (default star=1).
% G is the non-oriented graph and G1 the oriented result.
%
% see "Learning bayesian Networks from Data: A Efficient Approach Based on Information Theorie"
%     Jie Cheng, David Bell and Weird Liu.

% 0.
if nargin < 2, disp('Not enough arguments');return; end
if nargin < 3, node_sizes=max(data'); end
if nargin < 4, alpha = 0.05; end
if nargin < 5, star=1; end
G=G1;
N=length(G);
% reachability_matrix of G
M = expm(full(G)) - eye(length(G)); M = (M>0);

% 9.
fprintf('Thinning - separateA\n');

% Edges are examined in the inverse order of their Chi2 (or MI) score
s2(find(~G))=0;
[tmp ordre]=sort(s2(:));
ordre2=ordre(find(tmp>0));
[I J]=ind2sub([N N],ordre2);
ii=1:length(I);
for i=ii,
  %fprintf('%d-%d : ',I(i),J(i));
  G(I(i),J(i))=0;
  G(J(i),I(i))=0;
  trysep = try_to_separate_A(G,I(i),J(i),data,alpha,node_sizes);

  if ~trysep,
    G(I(i),J(i))=1;
    G(J(i),I(i))=1;
    %fprintf(' keep\n');
    %else
    %fprintf('delete\n');
  end
end

% 10.
fprintf('Thinning - separateB'); if star; fprintf('star'); end; fprintf('\n');
s2(find(~G))=0;
[tmp ordre]=sort(s2(:));
ordre2=ordre(find(tmp>0));
[I J]=ind2sub([N N],ordre2);
ii=1:length(I);
for i=ii,
  %fprintf('%d-%d : ',I(i),J(i));
  G(I(i),J(i))=0; G(J(i),I(i))=0;
  if star==0
    trysep = try_to_separate_B(G,I(i),J(i),data,node_sizes,alpha);
  else
    trysep = try_to_separate_B_star(G,I(i),J(i),data,node_sizes,alpha);
  end
  if ~trysep,
    G(I(i),J(i))=1; G(J(i),I(i))=1;
    %fprintf(' keep\n');
    %else
    %fprintf('delete\n');
  end
end

%11.
fprintf('Thinning - orient_edges\n');
U=G;
G=orient_edges(U,data,node_sizes,alpha);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I, N1, N2] = try_to_separate_A(G,node1,node2,data,alpha,node_sizes)
% [I, N1, N2] = try_to_separate_A(G,node1,node2,data,alpha,node_sizess)
%
% G is the current partially directed graph.
% I is a boolean : I=1 <==> separated.
% N1 : neighbors of node1 that are on an adjacency path between node1 and node2 (ditto for N2).
%
% see "Learning bayesian Networks from Data: A Efficient Approach Based on Information Theorie"
%     Jie Cheng, David Bell and Weird Liu.

% 0.
if node1==node2 | length(G)<2
  disp('Error: Check your arguments in try_to_separate_A.');I=0; return
end
N=size(data,1);
if nargin==4, alpha=0.05; node_sizes=max(data'); end
if nargin==5, node_sizes=max(data'); end

% 1.
N1=find(G(node1,:)==1);N01=N1;
GG1=G(setdiff(1:N,node1),setdiff(1:N,node1));
node22=node2-(node2>node1);
% reachability_matrix of GG1
M = expm(full(GG1)) - eye(length(GG1)); M = (M>0);
% N1 is the neighbors of node1 that are on the adjacency between node1 and node2
for i=N1
  j=i-(i>node1);
  if M(j,node22)~=1, N01=setdiff(N01,i); N1=setdiff(N1,i); end
  % 2.
  if ~G(node1,i), N1=setdiff(N1,i); end
end

N2=find(G(node2,:)==1);N02=N2;
GG2=G(setdiff(1:N,node2),setdiff(1:N,node2));
node12=node1-(node2<node1);
% reachability_matrix of GG2
M = expm(full(GG2)) - eye(length(GG2)); M = (M>0);
% N2 is the neighbors of node2 that are on the adjacency between node1 and node2
for i=N2
  j=i-(i>node2);
  if M(node12,j)~=1, N02=setdiff(N02,i); N2=setdiff(N2,i); end
  % 2.
  if ~G(node2,i), N2=setdiff(N2,i); end
end

%fprintf('%d : N1=',node1); fprintf('%d',N1); fprintf('\n');
%fprintf('%d : N2=',node2); fprintf('%d',N2); fprintf('\n');
% 3.
if length(N1)>length(N2), tmp=N1; N1=N2; N2=tmp; clear tmp, end
% 4.
C=N1;
for test=1:2
  if test==2, C=N2; end
  % 5.
  [I v1] = cond_indep_chisquare(node1,node2,C,data,'LRT',alpha,node_sizes);
  if I, %fprintf('%d-%d separated (%2.5f) by C=',node1,node2,v1); fprintf('%d',C); fprintf('\n');
    return,
    %else
    %fprintf('%d-%d not separated (%2.5f) by C=',node1,node2,v1); fprintf('%d',C); fprintf('\n');
  end;


  % 6.
  step6=1;
  while step6
    step6=0;
    if length(C)>=1
      v=[];
      for i=C
	Ci = setdiff(C,i);
	[I(i) v(i)] = cond_indep_chisquare(node1,node2,Ci,data,'LRT',alpha,node_sizes);
      end
      [vm ind] = min(v);

      % 7.
      if I(ind)
	I = 1; return
      else
	if vm < v1
	  v1 = vm;
	  C = setdiff(C,ind);
	  % goto step 6.
	  step6 = 1;
	end
      end
    end

    % 8.
    if test==2, I=0; return, end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = try_to_separate_B(G,node1,node2,data,node_sizes,alpha,N1,N2)
% I = try_to_separate_B(G,node1,node2,data,node_sizes,epsilon,N1,N2)
%
% G is the current partially directed graph.
% data(i,m), node i in case m.
% node_sizes is the vector of size of the attributs in data ( default max(data') )
% N1 (optionnal) is the neighbors of node1 that are on an adjacency path between node1 and node2 (ditto for N2).
% I is a boolean : I+1 <==> separated.
%
% see "Learning bayesian Networks from Data: A Efficient Approach Based on Information Theorie"
%     Jie Cheng, David Bell and Weird Liu.

% 0.
if node1==node2 | length(G)<2
  disp('Error: Verify your arguments in try_to_separate_B.');I=0; return
end
if nargin < 4,
  disp('Error : not enougth arguments'); I=0; return;
end
N=size(data,1);

% 1.
if nargin < 8
  N1=find(G(node1,:)==1);
  GG1=G(setdiff(1:N,node1),setdiff(1:N,node1));
  node22=node2-(node2>node1);
  M = expm(full(GG1)) - eye(length(GG1)); M = (M>0);
  for i=N1
    j=i-(i>node1);
    if M(j,node22)~=1, N1=setdiff(N1,i); end
  end
  N2=find(G(node2,:)==1);
  GG2=G(setdiff(1:N,node2),setdiff(1:N,node2));
  node12=node1-(node2<node1);
  M = expm(full(GG2)) - eye(length(GG2)); M = (M>0);
  for i=N2
    j=i-(i>node2);
    if M(node12,j)~=1, N2=setdiff(N2,i); end
  end
end
if nargin < 6, alpha=0.05; end
if nargin < 5, node_sizes=max(data'); end
M = expm(full(G)) - eye(length(G)); M = (M>0);

% 2.
N1b=[];
for i=N1
  NN1=find(G(i,:)==1);
  for j=NN1
    if M(i,j)~=1 & ~ismember(j,N1), N1b=union(N1b,j); end
  end
end

% 3.
N2b=[];
for i=N2
  NN2=find(G(i,:)==1);
  for j=NN2
    if M(i,j)~=1 & ~ismember(j,N2), N2b=union(N2b,j); end
  end
end

% 4.
if length(union(N1,N1b)) < length(union(N2,N2b))
  C=union(N1,N1b);
else
  C=union(N2,N2b);
end

% 5.
continu=1;
while continu
  l=length(C);
  %fprintf('%d',continu);
  [I v] = cond_indep_chisquare(node1,node2,C,data,'LRT',alpha,node_sizes);
  if I==1; return, elseif l<2, I=0; return, end

  % 6.
  Cb=C;
  for i=1:l
    Ci=setdiff(C,C(i));
    [I vi] = cond_indep_chisquare(node1,node2,Ci,data,'LRT',alpha,node_sizes);
    e = (v+1)/3;    % e is a small value...
	if I==1,return, elseif vi<v+e, Cb=setdiff(Cb,C(i)); end
  end

  % 7.
  if length(Cb) < l, C=Cb; else continu==0; I=0; return; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = try_to_separate_B_star(G,node1,node2,data,node_sizes,alpha,N1,N2)
% I = try_to_separate_B_star(G,node1,node2,data,epsilon,node_sizess)
%
% G is the current partially directed graph.
% I is a boolean.
% N1 is the neighbors of node1 that are on an adjacency path between node1 and node2 (ditto for N2).
%
% see "Learning bayesian Networks from Data: A Efficient Approach Based on Information Theorie"
%     Jie Cheng, David Bell and Weird Liu.

% 0.
if node1==node2 | length(G)<2
  disp('Error: Verify your arguments in try_to_separate_B_star.');I=0; return
end
if nargin < 4,
  disp('Error : not enougth arguments');I=0; return;
end
N=size(data,1);

% 1.
if nargin < 8
  N1=find(G(node1,:)==1);
  GG1=G(setdiff(1:N,node1),setdiff(1:N,node1));
  node22=node2-(node2>node1);
  M = expm(full(GG1)) - eye(length(GG1)); M = (M>0);
  for i=N1
    j=i-(i>node1);
    if M(j,node22)~=1, N1=setdiff(N1,i); end
    % 2.
    if ~G(node1,i), N1=setdiff(N1,i); end
  end
  N2=find(G(node2,:)==1);
  GG2=G(setdiff(1:N,node2),setdiff(1:N,node2));
  node12=node1-(node2<node1);
  M = expm(full(GG2)) - eye(length(GG2)); M = (M>0);
  for i=N2
    j=i-(i>node2);
    if M(node12,j)~=1, N2=setdiff(N2,i); end
    % 2.
    if ~G(node2,i), N2=setdiff(N2,i); end
  end
end
if nargin < 6, alpha=0.05; end
if nargin < 5, node_sizes=max(data'); end

% 3.
if length(N1)>length(N2), tmp=N1; N1=N2; N2=tmp; end

% 4.
C=N1;
l=length(C);
I=0;

% 5.
for test=1:2
  continu=1;
  %test
  if test==2 & ~isempty(N2)
    C=N2; l=length(C); IsInCi=zeros(1,l); IsInCi(l)=1;
  else IsInCi=zeros(1,l);
  end
  s=ones(1,l);
  % Pour tous les sous-ensemble Ci de C :
  while continu & ~isempty(C)
    Ci = setdiff(C.*IsInCi,0);
    [I vi] = cond_indep_chisquare(node1,node2,Ci,data,'LRT',alpha,node_sizes);
    if I, %fprintf('%d-%d separated (%2.5f) by C=',node1,node2,vi); fprintf('%d',Ci); fprintf('\n');
      return,
      %else
      %fprintf('%d-%d not separated (%2.5f) by C=',node1,node2,vi); fprintf('%d',Ci); fprintf('\n');
    end;
    %    if I, return, end

    if IsInCi==s, continu=0;
    else
      IsInCi(l)=IsInCi(l)+1;

      notOK=1; i=l;
      while notOK & i>1
	if IsInCi(i)>s(i), IsInCi(i)=0; IsInCi(i-1)=IsInCi(i-1)+1; else notOK=0; end
	i=i-1;
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G1 = orient_edges(G,data,node_sizes,alpha)
% G1 = orient_edges(G,data,node_sizes)
%
% G is the partially directed graph.
% data(i,m), node i in case m.
% node_sizes is the vector of size of the attributs in data ( default max(data') )
%
% see "Learning bayesian Networks from Data: A Efficient Approach Based on Information Theorie"
%     Jie Cheng, David Bell and Weird Liu.

% 0.
if nargin < 4, alpha=0.05; end
if nargin < 3, node_sizes=max(data'); end
if nargin < 2, disp(' Require at least two arguments.'); return, end
N=length(G);
G1=G;

% 1.
[Lnode1 Lnode2]=find(triu(1-triu(G),1)); %[Lnode1 Lnode2]=ind2sub([N N],find(triu(1-triu(G),1)));
for ii = 1:length(Lnode1),
  node1=Lnode1(ii);
  node2=Lnode2(ii);
  %fprintf('%d %d\n',node1,node2);
  N1 = find(G(node1,:)==1);
  N2 = find(G(node2,:)==1);
  if ~isempty(intersect(N1,N2))
    %fprintf('V1='); fprintf('%d ',N1); fprintf('\n');
    %fprintf('V2='); fprintf('%d ',N2); fprintf('\n');
    GG1 = G(setdiff(1:N,node1),setdiff(1:N,node1));
    node22 = node2-(node2>node1);
    % reachability_matrix of GG1
    M = expm(full(GG1)) - eye(length(GG1)); M = (M>0);
    % N1 is the neighbors of node1 that are on the adjacency between node1 and node2
    for i = N1
      j = i-(i>node1);
      if M(j,node22)~=1, N1=setdiff(N1,i);
      end
    end
    GG2 = G(setdiff(1:N,node2),setdiff(1:N,node2));
    node12 = node1-(node2<node1);
    % reachability_matrix of GG2
    M = expm(full(GG2)) - eye(length(GG2)); M = (M>0);
    % N2 is the neighbors of node2 that are on the adjacency between node1 and node2
    for i=N2
      j = i-(i>node2);
      if M(node12,j)~=1, N2=setdiff(N2,i); end
    end
    %fprintf('%d %d\n',node1,node2);
    %fprintf('N1='); fprintf('%d ',N1); fprintf('\n');
    %fprintf('N2='); fprintf('%d ',N2); fprintf('\n');

    % 2.
    M = expm(full(G)) - eye(length(G)); M = (M>0);
    N1b=N1;
    for i=N1
      NN1 = find(G(i,:)==1);
      for j=NN1
	if M(i,j)~=1 & ~ismember(j,N1), N1b=union(N1b,j); end
      end
    end
    %fprintf('N1''='); fprintf('%d ',N1b); fprintf('\n');

    % 3.
    N2b=N2;
    for i=N2
      NN2 = find(G(i,:)==1);
      for j=NN2
	if M(i,j)~=1 & ~ismember(j,N2), N2b=union(N2b,j); end
      end
    end
    %fprintf('N2''='); fprintf('%d ',N2b); fprintf('\n');

    % 4.
    if length(N1b) < length(N2b)
      C = N1b;
    else
      C = N2b;
    end
    %l=length(C);
    %fprintf('C='); fprintf('%d ',C); fprintf('\n');

    % 7.
    step5=1;
    while step5
      step5=0;
      %fprintf('.');
      % 5.
      l=length(C);
      [I v] = cond_indep_chisquare(node1,node2,C,data,'LRT',alpha,node_sizes);
      %fprintf('C='); fprintf('%d ',C);
      %fprintf(': %d %2.5f\n',I,v);

      step8=0;
      if I==1 & v~=0 % v < epsilon
	step8=1;
      else
	if  l==1
	  G1(C,node1)=0; G1(C,node2)=0;
	  fprintf('%d -> %d <- %d\n',node1,C,node2);
	  step8=1;
	end
      end
      %fprintf('%d\n',step8);

      % 6.
      if ~step8
	Cb=C;
	for i=1:l
	  Ci=setdiff(C,C(i));
	  [I vi] = cond_indep_chisquare(node1,node2,Ci,data,'LRT',alpha,node_sizes);
	  % e = (v+1)/3;    % e is a small value...
	  if I==1 % vi < v+e
	    Cb=setdiff(Cb,C(i));
	    if ismember(C(i),N1) & ismember(C(i),N2)
	      G1(C(i),node1)=0; G1(C(i),node2)=0;
	      fprintf('%d -> %d <- %d\n',node1,C(i),node2);
	    end
	    if I==1 % vi < epsilon
	      step8=1;
	    end
	  end
	end % for
      end % if

      % 7.
      if ~step8
	if length(Cb) < length(C), C=Cb; end
	if length(C) > 0, step5==1; end
      end
    end % while step5
    % step8 : passer � la paire de noeud suivant
    %fprintf('\n');
    %else
    %fprintf(' No common neighbor\n');
  end % if
end % for

% 11.
step9=0;
fprintf('Infering directions ');
test = pdag_to_dag(G1);
while ~isdag(test) %& step9<N
  step9=step9+1;
  %if ~isempty(test)
    %fprintf('.');
    % 9.
    for a=1:N, for b=1:N, for c=1:N,
	  if a~=b & b~=c & c~=a
	    if G1(a,b)==1 & G1(b,a)==0
	      %fprintf('%d -> %d ... \n',a,b);
	      if G1(b,c)==1 & G1(c,b)==1
		if G1(a,c)+G1(c,a)==0,
		  G1(c,b)=0;
		  fprintf('%d -> %d (9)\n',b,c);
		end
	      end
	    end
	  end
    end, end, end

    % 10.
    for a=1:N-1, for b=a+1:N
	if G1(a,b)==1 & G1(b,a)==1
	  GGG1=xor(G1,G1'); % matrice des arcs orient�s de G1
	  M = expm(full(GGG1)) - eye(length(GGG1)); M = (M>0);
	  if M(a,b)==1,
	    G1(b,a)=0;
	    fprintf('%d -> %d (10)\n',a,b);
	  end
	end
    end, end

      test = pdag_to_dag(G1);
      if isempty(test), G1=return_one_edge(G1); end
  %else
  %  G1, return,
  %end
end % while 11.
G1 = pdag_to_dag(G1);
fprintf('%d boucles\n',step9);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = isdag(G)
b = sum(sum(G.*G'));        % How many undirected arcs ? (x2)
b=~b & ~isempty(G);
if b
  M = expm(full(G)) - eye(length(G)); M = (M>0);
  b = b & find(sum(sum(eye(length(G)).*M))); % is there no cycle ?
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G=return_one_edge(G)
N=length(G); fam=[]; node2=[];
permnode = randperm(5); i=0; %fprintf('rev\n');
while isempty(fam) | isempty(node2) | i==N
 i=i+1;
 node = permnode(i); %node = ceil(rand(1)*N);
 fam=find(G(node,:)==1);
 par=find(G(:,node)==1);
 fam = myunion(fam, par);
 if ~isempty(fam)
  node2 = fam(ceil(rand(1)*length(fam)));
  par2=find(G(:,node)==1);
  if ~isempty(intersect(par2, node)), node2=[]; end
 end
end
if isempty(myintersect(node, par2)),
 G(node, node2)=0;
 G(node2, node)=1;
 fprintf('%d -> %d (Rev)\n',node, node2);
else
 G(node, node2)=1;
 G(node2, node)=0;
 fprintf('%d -> %d (Rev)\n',node2, node);
end
