function completed_pdag = complete_pattern(pdag)

%
% completed_dag = complete_pattern(pdag)
%
% uses Rules R1-R4 of Meek (1995) to complete
% orientations in a pdag as far as possible, 
% i.e. every compelled edge is oriented.
%
% (Rules R1-R4 are also summarized in Pearl (2000), p.51
%  and Neapolitan (2004), p. 546.)
%
% Since the PC algorithm also uses Rules R1-R3, their implementation was 
% copied (with some modifications) from function learn_struct_pdag_pc.
%
% Rule R4 is necessary here, since the orientations in the input 
% pdag do not just represent v-structures. 
%
% Imme Ebert-Uphoff (ebert@tree.com), 2007
%
		  
   DIAGNOSTICS_ON = false;

   n = length(pdag);
   old_pdag = zeros(n);
   %iter = 0;
   while ~isequal(pdag, old_pdag)
     %iter = iter + 1;
     old_pdag = pdag;

     % Rule R1
     [A,B] = find(pdag==-1); % a -> b
     for i=1:length(A)
       a = A(i); b = B(i);
       undirected = abs(pdag) + abs(pdag)';
       % Adjacency test in undirected matrix:
       %   a adjacent b  <=>  undirected(a,b) ==0
       % That's easier to use than adjacency test in pdag:
       %   a adjacent b  <=>  pdag(a,b)==0 and pdag(b,a)==0

       % Find all nodes c such that  b-c  and c not adjacent a
       C = find(pdag(b,:)==1 & undirected(a,:)==0); 
       if ~isempty(C)
         pdag(b,C) = -1; pdag(C,b) = 0; 
         if DIAGNOSTICS_ON
	     for j=1:length(C)   
	        fprintf('Rule 1: %d -> %d\n', b, C(j));
             end
          end
       end
     end

     % Rule R2
     [A,B] = find(pdag==1); % unoriented a-b edge
     for i=1:length(A)
       a = A(i); b = B(i);
       if any( (pdag(a,:)==-1) & (pdag(:,b)==-1)' ); 
         pdag(a,b) = -1; pdag(b,a) = 0; 
         if DIAGNOSTICS_ON
            fprintf('Rule 2: %d -> %d\n', a, b);
         end
       end
     end

     % Rule R3
     [A,B] = find(pdag==1); % a-b
     for i=1:length(A)
       a = A(i); b = B(i);
       C = find( (pdag(a,:)==1) & (pdag(:,b)==-1)' );
       % C contains nodes c s.t. a-c->b-a

       % Extract lines and columns corresponding only to the set of nodes C
       core = pdag(C,C);

       % Prepare adjacency test:
       unoriented = abs(core) + abs(core)';  
       % Now:  a non-adjacent b <==> unoriented(a,b) == 0

       % Prepare to detect existence of non-adjacent pairs of nodes in C.
       % Set diagonal to 1, to prevent finding pairs of IDENTICAL nodes:
       unoriented = setdiag(unoriented, 1);
       if any(unoriented(:)==0) % C contains 2 different non adjacent elements
         pdag(a,b) = -1; pdag(b,a) = 0; 
         if DIAGNOSTICS_ON
            fprintf('Rule 3: %d -> %d\n', a, b);
         end
       end
     end

     % Rule 4
     [A,B] = find(pdag==1); % unoriented a-b edge
     for i=1:length(A)
       a = A(i); b = B(i);

       % Prepare adjacency test:
       % unoriented(i,j) is 0 (non-adj) or 1 (directed) or 2 (undirected)
       unoriented = abs(pdag) + abs(pdag)';

       % Find c such that c -> b and a,c are adjacent (a-c or a->c or a<-c) 
       C = find( (pdag(:,b)==-1)' & (unoriented(a,:)>=1) );  
       for j=1:length(C)
          c = C(j);
          % Check whether there is any node d, such that
          % d->c  AND  a-d  AND  b NOT adjacent to d
          if any( (pdag(:,c)==-1)' & (pdag(a,:)==1) & (unoriented(b,:)==0) )
	     pdag(a,b) = -1;  pdag(b,a) = 0;  
             if DIAGNOSTICS_ON
                fprintf('Rule 4: %d -> %d\n', a, b);
             end
          end
       end
     end

   end % end of while

   % Oriented all possible edges.  Return result.
   completed_pdag = pdag;
end
  
