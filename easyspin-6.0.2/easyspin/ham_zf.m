% ham_zf  Electronic zero field interaction Hamiltonian 
%
%   F = ham_zf(SpinSystem)
%   F = ham_zf(SpinSystem,Electrons)
%   F = ham_zf(SpinSystem,Electrons,'sparse')
%
%   Returns the electronic zero-field interaction (ZFI)
%   Hamiltonian of the system SpinSystem, in units of MHz.
%
%   If the vector Electrons is given, the ZFI of only the
%   specified electrons is returned (1 is the first, 2 the
%   second, etc). Otherwise, all electrons are included.
%
%   If 'sparse' is given, the matrix is returned in sparse format.
