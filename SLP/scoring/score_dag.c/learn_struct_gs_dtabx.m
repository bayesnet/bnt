function [dag,best_score] =	learn_struct_gs_dtabx(data, nodesizes, seeddag, varargin)
%
% LEARN_STRUCT_GS(data,seeddag)	learns a structure of Bayesian net by Greedy Search.
% dag =	learn_struct_gs(data, nodesizes, seeddag)
%
% dag: the final structure matrix
% Data : training data,	data(i,m) is the m obsevation of node i
% Nodesizes: the size array	of different nodes
% seeddag: given seed Dag for hill climbing, optional
%
% by Gang Li @ Deakin University (gli73@hotmail.com)
%
% -----------------------------------------------------
%
% Modified from	learn_struct_gs (SLP 1.3)
% to learn structure of	BN with	tabular	nodes:
%
% 1) make use of score decomposibility 
% 2) replace score_family with score_family_x.c that calculates 
%    Bayesian score with default parameters (non-adjustable)
%
% by Darima	<darrimma@yahoo.com>, 28/12/2005
%

[N ncases] = size(data);
if (nargin < 3 ) 
	seeddag	= zeros(N,N); %	mk_rnd_dag(N); %call BNT function
elseif ~acyclic(seeddag)
	seeddag	= mk_rnd_dag(N); %zeros(N,N);
end;

% set default params (the same as in score_dags)
for i=1:N
  type{i} = 'tabular';
  params{i} = { 'prior_type', 'dirichlet', 'dirichlet_weight', 1 };
end
scoring_fn = 'bayesian';
discrete = 1:N;
verbose	 = 'yes';

% get params
args = varargin;
nargs =	length(args);
if length(args)	> 0
	if isstr(args{1})
		for	i =	1:2:nargs
			switch args{i}
			case 'scoring_fn', scoring_fn =	args{i+1};
            case 'type',       type = args{i+1}; 
            case 'discrete',   discrete = args{i+1}; 
            case 'params',     
                if isempty(args{i+1}), params = cell(1,n); 
                else params = args{i+1};  end;           
			case 'verbose',	 verbose  =	strcmp(args{i+1},'yes');
			end;
		end;
	end;
end;

%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
datax = data-ones(N,ncases);     %%+++++++++++++++++++++++ IMPORTANT!!!!
%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

best_score = score_dags(data,nodesizes,	{seeddag},'scoring_fn',scoring_fn);
done = 0;
it = 0;

while ~done
	[dags,op,nodes]	= mk_nbrs_of_dag(seeddag);
	nbrs = length(dags);
%------ DEBUG    
    fprintf('DEBUG: nbrs = %6d\n',nbrs);
    time_inloop = cputime;
%------ DEND
%     scores	= score_dags(data, nodesizes, dags,'scoring_fn',scoring_fn);
	scores = zeros(nbrs,1);
	for	i =	1:nbrs
		xj = nodes(i,2);	
		ps_old = parents(seeddag, xj)';
		ps_new = parents(dags{i}, xj)';
        scor_old = score_family_x([datax(ps_old,:);datax(xj,:)],...
            [nodesizes(ps_old),nodesizes(xj)]);
        scor_new = score_family_x([datax(ps_new,:);datax(xj,:)],...
            [nodesizes(ps_new),nodesizes(xj)]);
		scores(i) =	best_score - scor_old +	scor_new;
		if isequal(op{i},'rev')
		    xi = nodes(i,1);
			ps_old = parents(seeddag, xi);
			ps_new = parents(dags{i}, xi);
            scor_old = score_family_x([datax(ps_old,:);datax(xi,:)],...
                [nodesizes(ps_old),nodesizes(xi)]);
            scor_new = score_family_x([datax(ps_new,:);datax(xi,:)],...
                [nodesizes(ps_new),nodesizes(xi)]);
			scores(i) =	scores(i) - scor_old +	scor_new;            
		end
    end    
	max_score =	max(scores);
	new	= find(scores == max_score );
%------ DEBUG
    fprintf('       -> max_score = %7.5f\n',max_score);
    fprintf('       -> find(scores == max_score): %d... of %d\n',...
        new(1),length(new));
%------ DEND
	if ~isempty(new) & (max_score >	best_score)
		p =	sample_discrete(normalise(ones(1, length(new))));
		best_score = max_score;
		seeddag	= dags{new(p)};
	else
		done = 1;
	end;    
    it = it+1;
%------ DEBUG
    time_inloop = cputime-time_inloop;
    fprintf('       time = %12.5f\n',time_inloop);
%------ DEND
end;
dag	= seeddag;

%----------------------------------------------------------------------

outcount = 0; 
best_score = score_dags(data,nodesizes,	{seeddag},'scoring_fn',scoring_fn);
while outcount < 2
	innercount = 0;
	for	i=1:N
		for	j=1:N
		   if i==j,	continue;	 end;
		   if seeddag(i,j) == 0	 % No edge i-->j, then try to add it
			   tempdag = seeddag;
			   tempdag(i,j)	= 1;
			   if acyclic(tempdag)
					temp_score = score_dags(data,nodesizes,	{tempdag},'scoring_fn',scoring_fn);
					if temp_score >	best_score
						seeddag	= tempdag;
						best_score=	temp_score;
						innercount = innercount	+1;
					end;
			   end
		   else	 % exists edge i--j, then try reverse it or	remove it
			   tempdag = seeddag;
			   tempdag(i,j)	= 0; tempdag(j,i) =	1; 
			   if acyclic(tempdag)
				   temp_score =	score_dags(data,nodesizes, {tempdag},'scoring_fn',scoring_fn);
				   if temp_score > best_score
					   seeddag = tempdag;
					   best_score =	temp_score;
					   innercount =	innercount +1;
				   else
					   tempdag = seeddag;
					   tempdag(i,j)	= 0;
					   temp_score =	score_dags(data,nodesizes, {tempdag},'scoring_fn',scoring_fn);
					   if temp_score > best_score
						   seeddag = tempdag;
						   best_score= temp_score;
						   innercount =	innercount +1;
					   end;
				   end;
			   else
				   tempdag = seeddag;
				   tempdag(i,j)=0;
				   temp_score =	score_dags(data,nodesizes, {tempdag},'scoring_fn',scoring_fn);
				   if temp_score > best_score
					   seeddag = tempdag;
					   best_score= temp_score;
					   innercount =	innercount +1;
				   end;
			   end;
		   end;
		end; % end for j
	end; % end for i
	if innercount == 0
		outcount = outcount	+1;
	end;
end;  %	end	while

%------ DEBUG
fprintf('DEBUG: Number of iterations = %d\n',it);
fprintf('DEBUG: Outcount = %d\n',outcount);
%------ DEND
dag	= seeddag;

