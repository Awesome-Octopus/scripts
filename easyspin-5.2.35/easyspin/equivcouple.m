% equivcouple   Combination of equivalent spins 
%
%  [F,N] = equivcouple(I,n)
%
%  The states due to n spins-I can be combined to give
%  a set of independent spins. Their quantum numbers
%  are returned in F, their respective abundances in N.
%
%  Example:
%
%     [F,N] = equivcouple(1/2,5)
%
%  5 spins-1/2 give rise to a first-order splitting pattern
%  [1 5 10 10 5 1] (see the function equivsplit). This can be
%  decomposed into one spin-5/2, four spin-3/2 and five spin-1/2
%  according to
%
%          5  5          5 spins-1/2
%       4  4  4  4       4 spins-3/2
%    1  1  1  1  1  1    1 spin-5/2
%   ------------------
%    1  5 10 10  5  1    sum
%
%  so F = [2.5 1.5 0.5] and N = [1 4 5].
%
%  In group theoretical terms, this corresponds to the reduction
%  of a product of n irreps of dimension 2*I of the rotation group
%  into a direct sum of irreps (Clebsch-Gordan decomposition).
