% soint  Spin-orbit interaction Hamiltonian 
%
%   H = soint(SpinSystem)
%   H = soint(SpinSystem,Spins)
%   H = soint(SpinSystem,Spins,'sparse')
%
%   Returns the spin orbit interaction (SOI)
%   Hamiltonian, in MHz.
%
%   Input:
%   - SpinSystem: Spin system structure. SOI
%       parameters are the spin-orbit coupling in soc 
%       and the orbital reduction facotr in orf.
%   - Spins: If given, specifies spins
%       for which the SOI should be computed. If
%       absent, all electrons are included.
%   - 'sparse': If given, the matrix is returned in sparse format.
