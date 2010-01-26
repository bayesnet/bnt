function pot = recursive_combine_pots(pot1, pot2)
% RECURSIVE_COMBINE_POTS recursive combine two potentials
% pot = recursive_combine_pots(pot1, pot2)

pot1 = reduce_pot(pot1);
pot2 = reduce_pot(pot2);
% Recursion is stopped, if recusive-combination is defined by direct combination, 
% i.e. if the domain of one potential is disjoint from the head of the other.
if (isempty(myintersect(pot1.domain,pot2.cheaddom))|...
    isempty(myintersect(pot1.cheaddom,pot2.domain)))    
    pot = direct_combine_pots(pot1,pot2);
else 
    % Test wether one of the set-differences is not empty 
    % as defined in Lauritzen99 "Stable Local Computation with Conditional Gaussian Distributions"
    % on page 9
    D12 = mysetdiff(pot1.cheaddom, pot2.domain);
    D21 = mysetdiff(pot2.cheaddom, pot1.domain);
    if (isempty(D12) & isempty(D21))
       assert(0,'Recursive combination is not defined');
    end

    if ~isempty(D12)
        % Calculate the complementary potential for the set 
        % D1\D12 as defined in Lauritzen 99, page 9
    keep = mysetdiff(pot1.domain,D12);
        [margpot, comppot] = complement_pot(pot1,keep);
        margpot = reduce_pot(margpot);
        comppot = reduce_pot(comppot);
        pot = direct_combine_pots( recursive_combine_pots(margpot, pot2), comppot);
    elseif ~isempty(D21)
        keep = mysetdiff(pot2.domain,D21);
        [margpot, comppot] = complement_pot(pot2,D21);
        margpot = reduce_pot(margpot);
        comppot = reduce_pot(comppot);
        pot = direct_combine_pots( recursive_combine_pots(pot1, margpot), comppot);
    end
end



