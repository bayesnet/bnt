function mult_by_table_global(bigT, bigdom, bigsz, smallT, smalldom, smallsz)

% all arguments are read only
global NEWBIGT_GLOBAL

Ts = extend_domain_table(smallT, smalldom, smallsz, bigdom, bigsz);
NEWBIGT_GLOBAL = bigT(:) .* Ts(:);
