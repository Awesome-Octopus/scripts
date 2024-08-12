% oripotentialplot     Plot orientational potential
%
%    oripotentialplot(Potential)
%
%   Plots the orientational potential U(alpha,beta,gamma) specified in Potential.
%
%   Input:
%     Potential    array with L, M, K, and lambda values for the Wigner
%                  expansion representation of the potential:
%                  [L1 M1 K1 lambda; L2 M2 K2 lambda2; ...]
%
%   Example:
%     Potential = [2 2 1 1.4]; % L=2, M=2, K=1, lambda=1.4
%     oripotentialplot(Potential);
