% stev  Extended Stevens spin operator matrices
%
%   Op = stev(S,[k,q])
%   Op = stev(S,[k,q,iSpin])
%   Op = stev(__,'sparse')
%
%   Constructs extended Stevens operator matrices of rank k and component q
%   for spin S. Rank and component must satisfy 0<=k<=2*S and -k<=q<=k.
%
%   k is the rank of the operator, and q is the component. q<0 correspond to
%   the sin tesseral components O_k^q(s), and q>0 to the cos tesseral components
%   O_k^q(c). Ranks from 0 to 12 are supported. The most common ones are 2, 4 and 6.
%
%   If S is a vector representing the spins of a spin system, Op is computed for
%   the spin number iSpin (e.g. the second if iSpin==2) in the state space of
%   the full spin system.
%
%   The extended Stevens operators are tesseral (as opposed to spherical) tensor
%   operators and are all Hermitian.
%
%   Input:
%   - S: spin quantum number, or vector thereof
%   - k,q: indices specifying O_k^|q|(c) for q>=0 and O_k^|q|(s) for q<0
%   - iSpin: index of the spin in the spin vector for
%       which the operator matrix should be computed
%   - 'sparse': if given, return matrix in sparse and not in full format
%
%   Output:
%   - Op: extended Stevens operator matrix
%
%   Examples:
%    Get O_4^2(c) for a spin 5/2:
%       stev(5/2,[4,2])
%    Get O_6^5(c) for the second spin in a spin system with two spins-3/2:
%       stev([3/2 3/2],[6,5,2])
