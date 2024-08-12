% equivsplit   Equivalent nuclei: EPR splitting pattern 
%
%  Ampl = equivsplit(I,n)
%
%  Computes the line intensity pattern of an EPR
%  spectrum due to an S=1/2 and n equivalent nuclear spins
%  with quantum number I. First-order perturbations
%  are assumed.
%
%  Example:
%
%     Intensity = equivsplit(1/2,5)
%
%  5 spins-1/2 give rise to a first-order splitting
%  pattern Intensity = [1 5 10 10 5 1] according to
%
%              1              0 spins
%            1   1            1 spin-1/2
%          1   2   1          2 spins-1/2
%        1   3   3   1        3 spins-1/2
%      1   4   6   4   1      4 spins-1/2
%    1   5  10  10   5   1    5 spins-1/2
