%  data generation
clear all
close all

ddd = datestr(now);
ddd([12 15 18])='-' ;
fnd=[ddd '.txt'];
diary(fnd)

% Metaparams
hyp = 'mar';
ntests = 2;   % change it !
nessais = 2;  % change it !
m = 1000;     % change it !

% INIT
missing = 2;
misv = -9999;
%test = zeros(nessais,6);
test = zeros(nessais,ntests);
taille = zeros(nessais,ntests);
kiki = zeros(nessais,1);

%%%%%%%%%% TESTS %%%%%%%%%%%
for essai=1:nessais

  N = ceil(rand*11) + 2
  base_proba = .15 + rand/4;

  fan_in = ceil(rand*4)
  dag=mk_rnd_dag(N,fan_in)
  discrete=1:N;
  for i=1:N, node_sizes(i) = ceil(rand*4)+1; end
  bnet = mk_bnet(dag, node_sizes, discrete);
  node_sizes
  for i=1:N, bnet.CPD{i} = tabular_CPD(bnet, i); end

  if prod([hyp(1:3) == 'mca']+0),
    bnet_miss = gener_MCAR_net(bnet, base_proba);
    disp(' - MCAR net generated');
  elseif prod([hyp(1:3) == 'mar']+0),
    bnet_miss = gener_MAR_net(bnet, base_proba);
    disp(' - MAR net generated');
  end

  [data, comp_data, bnet_miss, taux, bnet_orig, notok, dT] = gener_data_from_bnet_miss(bnet_miss, m, base_proba ,1);
  base_proba

  kiki(essai) = 1-chi2cdf(dT,1);

  %test(essai,1)=1-chi2cdf(d,df);
  node0=[]; param0=[];
  for i=1:ntests
    % choose a node
    node=ceil(rand*N);
    while ~isempty(myintersect(node,node0)) & i<N-2,
      node = ceil(rand*N);
    end
    node0 = myunion(node0,node);

    % choose a param
    par = parents(bnet_miss.dag, node+N);
    sorted=[]; ord=[];
    for par1=1:length(par)
      ord(par1)=find(par(par1)==bnet_miss.order);
    end
    [tmp, sorted] = sort(ord);
    family = par(sorted);

    state = []; stateM = [];
    for par1=1:length(family)
      state(par1) = ceil(rand*bnet_miss.node_sizes(family(par1)));
      if state(par1) == bnet_miss.node_sizes(family(par1)), stateM(par1) = misv; else, stateM(par1) = state(par1); end
    end

    if prod([hyp(1:3) == 'mca']+0),
      family = family-N;
    elseif prod([hyp(1:3) == 'mar']+0),
      family = family-2*N;
    end
    family(end+1)=node;
    state(end+1) = missing; stateM(end+1)=misv;

    CPT = CPT_from_bnet(bnet_miss);
    value = CPT{node+N}(subv2ind(bnet_miss.node_sizes(family), state))

    % test this param "value"
    datamat = bnt_to_mat(data(family,:), misv);
    indj=1:m;
    for par1 = 1:(length(family)-1)
      [tmp, indj] = find(datamat(par1,indj)==stateM(par1));
    end
    mtest = length(indj);
    [tmp, indj] = find(datamat(length(family),indj)==stateM(length(family)));

    toto=mtest; nbr_miss=length(indj);
d = ((nbr_miss-value*toto)^2)/(value*toto)...
  + (((toto-nbr_miss)-(1-value)*toto)^2)/((1-value)*toto)

    family
    %df = prod(bnet_miss.node_sizes(family)-1);
    df = 1; % present or missing
    %for chiv=1:nchi,
    %  if d>chi2_table(chi2(chiv),df), test(essai,i,chiv) = 1; end
    %end
    cumchi=chi2cdf(d,df)
    test(essai,i) = 1-cumchi;
    taille(essai,i) = toto;

    %fprintf('%2.1f%% OF MISSING DATA (%2.1f%%, Khi2 : %2.1f > %2.1f)\n', round(base_proba*10000)/100, round(taux*10000)/100, d, chi2(nchi));

  end
  clear bnet_miss bnet datamat data node0 family node_sizes CPT indj tmp
end

kiki
taille
test

%eval(['save resgener' ddd 'kiki taille test hyp ntests nessais = 50 m'])

diary off

