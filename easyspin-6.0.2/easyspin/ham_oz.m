% ham_ez  Orbital Zeeman interaction Hamiltonian 
%
%   H = ham_oz(SpinSystem, B)
%   H = ham_oz(SpinSystem, B, oam)
%   H = ham_oz(SpinSystem, B, oam, 'sparse')
%   [mux, muy, muz] = ham_oz(SpinSystem)
%   [mux, muy, muz] = ham_oz(SpinSystem, oam)
%   [mux, muy, muz] = ham_oz(SpinSystem, oam, 'sparse')
%
%   Returns the orbital Zeeman interaction Hamiltonian matrix for
%   the orbital angular momenta 'oam' of the spin system 'SpinSystem'.
%
%   Input:
%   - SpinSystem: Spin system structure.
%   - B: Magnetic field vector, in millitesla.
%   - oam: Vector of indices for orbital angular momenta to include. If oam
%     is omitted, all orbital angular momenta are included.
%   - 'sparse': If given, results returned in sparse format.
%
%   Output:
%   - mux, muy, muz: Components of the magnetic dipole moment operator
%     for the selected orbital angular momenta as defined by mui=d(H)/d(B_i)
%     i=x,y,z where B_i are the cartesian components of
%     the external field in the molecular frame. Units are MHz/mT = GHz/T.
%     To get the full orbital Zeeman Hamiltonian, use
%     H = -(mux*B(1)+muy*B(2)+muz*B(3)), where B is the magnetic field vector in mT.
%   - H: Orbital Zeeman Hamiltonian matrix.
