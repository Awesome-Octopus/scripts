% zeemanho  Multiple Order Zeeman interaction Hamiltonian 
%
%
%   H = zeemanho(SpinSystem, B)
%   H = zeemanho(SpinSystem, B, Spins)
%   H = zeemanho(SpinSystem, B, Spins, 'sparse')
%   H = zeemanho(SpinSystem, B, Spins, 'sparse', lB)
%   H = zeemanho(SpinSystem, B, [], '', lB)
%
%
%   [G0,G1] = zeemanho(SpinSystem)
%   [G0,G1,G2] = zeemanho(SpinSystem)
%   [G0,G1,G2,G3] = zeemanho(SpinSystem)
%   GlB = zeemanho(SpinSystem, Spins, 'sparse',lB)
%   [G0,...] = zeemanho(SpinSystem, Spins, 'sparse')
%   cG = zeemanho(SpinSystem)
%   cG = zeemanho(SpinSystem, [], Spins, 'sparse',lB)
%   (for cell output and spin selection, empty field for B required)
%
%   Returns the multiple order Zeeman interaction Hamiltonian for
%   the electron spins 'Spins' of the spin system 'SpinSystem'.
%
%   Input:
%   - SpinSystem: Spin system structure.
%   - B: Magnetic field vector, in milliTesla. If B is ommited zero field
%   is used (for single output).
%   - Spins: Vector of spin numbers.  If Spins is omitted,
%     all electron spins are included.
%   - 'sparse': If given, results returned in sparse format.
%   - lB: If given only terms of order lB in the magnetic field are
%   returned
%
%   Output:
%   - H: the Hamiltonian of the Zeeman interaction.
%   -[G0,G1,G2,G3]: components of the multiple order Zeeman interaction Hamiltonian
%     for the selected spins, defined by the 0th, 1th, 2nd, and 3rd derivative of
%     of the parts of the Hamiltonian which contains the magnetic up to
%     this order.  Gn contain (d^n/dB^^n) H, a tensor of rank n. 
%     Units are MHz/(mT)^n.
%     So  G3{1,2,3} contain d^3 /(dB_x dB_y dB_z) H and is in units of 
%     MHz/(mT)^3 = 10^15 Hz/T^3. 
%   -cG cell containing the G
%
%   Uses the Hamiltonian as given in
%   McGavin, Tennant and Weil, 	J. Magn. Reson. 87,92-109 (1990)
%
%   it is complete in the sense that it contain all usual terms (in
%   principal also the nuclear and hyperfine for a single nuclei, but this
%   is not implemented) for a single spin, but in a non-popular convention.
