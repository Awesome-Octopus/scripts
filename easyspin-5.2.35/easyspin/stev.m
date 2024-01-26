% stev  Extended Stevens spin operator matrices
%
%   Op = stev(S,k,q)
%   Op = stev(S,k,q,iSpin)
%   Op = stev(S,k,q,iSpin,'sparse')
%
%   Constructs extended Stevens operator matrices for
%   0<=k<=2*S and -k<=q<=k for the spin S.
%
%   If S is a vector representing the spins of a
%   spin system, Op is computed for the spin number
%   iSpin (e.g. the second if iSpin==2) in the state
%   space of the full spin system. If iSpin is omitted,
%   iSpin is set to 1.
%
%   All k values from 0 to 12 are supported. The
%   most common ones are 2, 4 and 6.
%
%   The extended Stevens operators are tesseral (as
%   opposed to spherical) tensor operators and are
%   therefore all Hermitian.
%
%   Input:
%   - S: spin quantum number, or vector thereof
%   - k,q: indices specifying O_k^q
%   - iSpin: index of the spin in the spin vector for
%       which the operator matrix should be computed
%
%   Output:
%   - Op: extended Stevens operator matrix
%
%   Examples:
%    To obtain O_4^2 for a spin 5/2, type
%       stev(5/2,4,2)
%    To obtain O_6^5 for the second spin in a spin
%    system with two spins-3/2, type
%       stev([3/2 3/2],6,5)
