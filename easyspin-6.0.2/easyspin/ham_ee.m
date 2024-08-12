% ham_ee  Electron-electron spin interaction Hamiltonian 
%
%   F = ham_ee(SpinSystem)
%   F = ham_ee(SpinSystem,eSpins)
%   F = ham_ee(SpinSystem,eSpins,'sparse')
%
%   Returns the electron-electron spin interaction (EEI)
%   Hamiltonian, in MHz.
%
%   Input:
%   - SpinSystem: Spin system structure. EEI
%       parameters are in the ee, eeFrame, and ee2 fields.
%   - eSpins: If given, specifies electron spins
%       for which the EEI should be computed. If
%       absent, all electrons are included.
%   - 'sparse': If given, the matrix is returned in sparse format.
