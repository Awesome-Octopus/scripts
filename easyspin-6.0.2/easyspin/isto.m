% isto  Irreducible spherical tensor operators (ISTOs)
%
%   T = isto(J,kq)
%   T = isto(J,kqi)
%   T = isto(___,'sparse')
%
% Input:
%   J    array of angular momentum quantum numbers (0, 1/2, 1, etc.);
%   kq   array of ranks k and projections q for each of the angular momenta in J
%        each row contains [k q] for the corresponding J
%   kqi  array of ranks k, projections q, and ang.mom. indices i; each row
%        contains [k q i]. For any ang.mom. not listed, k=0 and q=0 are assumed.
%
% Output:
%   T    irreducible spherical tensor operator matrix
%
% If 'sparse' is given, the matrix is returned in sparse storage format,
% otherwise it is returned in full format.
%
% Examples:
%   T = isto(3/2,[2 0])            % T^2_0 for a spin 3/2
%   T = isto([1/2 1/2],[1 1 2])    % T^1_1 for the second spin
