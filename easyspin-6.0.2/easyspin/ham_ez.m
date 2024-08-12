% ham_ez  Electron Zeeman interaction Hamiltonian 
%
%   H = ham_ez(SpinSystem, B)
%   H = ham_ez(SpinSystem, B, eSpins)
%   H = ham_ez(SpinSystem, B, eSpins, 'sparse')
%   [mux, muy, muz] = ham_ez(SpinSystem)
%   [mux, muy, muz] = ham_ez(SpinSystem, eSpins)
%   [mux, muy, muz] = ham_ez(SpinSystem, eSpins, 'sparse')
%
%   Returns the electron Zeeman interaction Hamiltonian matrix for
%   the electron spins 'eSpins' of the spin system 'SpinSystem'.
%
%   Input:
%   - SpinSystem: Spin system structure.
%   - B: Magnetic field vector, in millitesla.
%   - eSpins: Vector of indices for electron spins to include. If eSpins is
%     omitted, all electron spins are included.
%   - 'sparse': If given, results returned in sparse format.
%
%   Output:
%   - mux, muy, muz: Components of the magnetic dipole moment operator
%     for the selected electron spins as defined by mui=-d(H)/d(B_i)
%     i=x,y,z where B_i are the cartesian components of
%     the external field in the molecular frame. Units are MHz/mT = GHz/T.
%     To get the full electron Zeeman Hamiltonian, use
%     H = -(mux*B(1)+muy*B(2)+muz*B(3)), where B is the magnetic field vector
%     in mT.
%   - H: Electron Zeeman Hamiltonian matrix.
