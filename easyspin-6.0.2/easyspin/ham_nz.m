% ham_nz  Nuclear Zeeman interaction Hamiltonian 
%
%   H = ham_nz(SpinSystem, B)
%   H = ham_nz(SpinSystem, B, nSpins)
%   H = ham_nz(SpinSystem, B, nSpins, 'sparse')
%   [mux, muy, muz] = ham_nz(SpinSystem)
%   [mux, muy, muz] = ham_nz(SpinSystem, nSpins)
%   [mux, muy, muz] = ham_nz(SpinSystem, nSpins, 'sparse')
%
%   Returns the nuclear Zeeman interaction Hamiltonian matrix for
%   the nuclear spins 'nSpins' of the spin system 'SpinSystem'.
%
%   Input:
%   - SpinSystem: Spin system structure.
%   - B: Magnetic field vector, in millitesla.
%   - nSpins: Vector of indices for nuclear spins to include. If nSpins is
%     omitted, all nuclear spins are included.
%   - 'sparse': If given, results returned in sparse format.
%
%   Output:
%   - mux, muy, muz: Components of the magnetic dipole moment operator
%     for the selected nuclear spins as defined by mui=-d(H)/d(B_i)
%     i=x,y,z where B_i are the cartesian components of
%     the external field in the molecular frame. Units are MHz/mT = GHz/T.
%     To get the full nuclear Zeeman Hamiltonian, use
%     H = -(mux*B(1)+muy*B(2)+muz*B(3)), where B is the magnetic field vector in mT.
%   - H: Nuclear Zeeman Hamiltonian matrix.
