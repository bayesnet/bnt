function bigT = mult_by_table2(bigT, bigdom, bigsz, smallT, smalldom, smallsz)
% MULT_BY_TABLE 
% bigT = mult_by_table(bigT, bigdom, bigsz, smallT, smalldom, smallsz)
%

%Ts = extend_domain_table(smallT, smalldom, smallsz, bigdom, bigsz);
%bigT(:) = bigT(:) .* Ts(:); % must have bigT(:) on LHS to preserve shape

% extend_domain_table has a lot of overhead for small tables,
% since it calls myreshape and myrepmat, which check for 1 dimensional case.
% Here, we check up front.

if length(bigdom)==1 % vector
  bigT = bigT .* smallT; % smallT can be scalar or vector
else
  if (length(bigsz) == length(smallsz)) & all(bigsz == smallsz)
    bigT = bigT .* smallT;
  else
    map = find_equiv_posns(smalldom, bigdom);
    sz = ones(1, length(bigdom));
    sz(map) = smallsz;
    smallT = reshape(smallT, sz); % add dimensions of size 1 for missing domain
    % we can use reshape instead of myreshape, because we know length(sz)>1
    sz = bigsz;
    sz(map) = 1; % don't replicate along small domain, which is shared
    % we can use repmat instead of myrepmat, because we know length(sz)>1
    smallT = repmat(smallT, sz(:)');
    bigT(:) = bigT(:) .* smallT(:);
  end
end
