function pot = combine_pots(pot1, pot2)
% COMBINE_POTS combine two potentials 
% pot = combine_pots(pot1, pot2)

% Reduce both potentials before trying to combine them. 
% Cf. "Stable Local computation with Conditional Gaussian Distributions", page 9
% Consider again two potentials with minimal tail

% Guarantee minimal tails. If pot1 or pot2 are minimal, they are not changed
pot1 = reduce_pot(pot1);
pot2 = reduce_pot(pot2);

%if the intersect set of these two potentials' head conts. combination is undifined
if ~isempty( myintersect(pot1.cheaddom, pot2.cheaddom) )
    return;
end

if  isempty( myintersect(pot1.domain, pot2.cheaddom) ) | isempty( myintersect(pot2.domain, pot1.cheaddom))
    % if satisfy the condition of directed combine
    pot = direct_combine_pots(pot1, pot2);
else
    % perform recursive combine
    pot = recursive_combine_pots(pot1, pot2);
end
