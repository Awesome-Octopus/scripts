% isot2stev  Transformation from ISTOs to Stevens operators
%   
%   C = isto2stev(k)
%
% Provide transformation matrix to transform a rank-k irreducible spherical tensor
% operator (ISTO) T_k to a vector of rank_k extended Steven operators O_k with
%
%   O_k = C *T_k
%
% T_k is the ISTO of rank k with 2*k+1 components T(k,q=k), T(k,k-1), etc.
% and O_k is the vector of Stevens operators O(k,q=k), O(k,k-1), ... O(k,0),...
% O(k,-k). Stevens operators with non-negative q are cosine tesseral operators,
% and those with negative q are sine tesseral operators.
%
% C is the (2k+1)x(2k+1) transformation matrix. Its only non-negative elements
% are on the diagonal and the antidiagonal.
%
% For the inverse transform, use inv(C): T_k = inv(C)*O_k
