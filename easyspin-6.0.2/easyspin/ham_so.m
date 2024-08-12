% ham_so  Spin-orbit interaction Hamiltonian 
%
%   H = ham_so(SpinSystem)
%   H = ham_so(SpinSystem,eSpins)
%   H = ham_so(SpinSystem,eSpins,'sparse')
%
%   Returns the spin orbit interaction (SOI) Hamiltonian, in MHz.
%
%   Input:
%   - SpinSystem: Spin system structure. The spin-orbit coupling is
%       in SpinSystem.soc.
%   - Spins: If given, specifies electron spins for which the SOI should be
%       computed. If absent, all electrons are included.
%   - 'sparse': If given, the matrix is returned in sparse format.
%
%   Output:
%   - H: Hamiltonian matrix containing the SOI for spins specified in eSpins.
