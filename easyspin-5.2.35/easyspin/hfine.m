% hfine  Hyperfine interaction Hamiltonian 
%
%   F = hfine(System)
%   F = hfine(System,Spins)
%   F = hfine(System,Spins,'sparse')
%
%   Returns the hyperfine interaction Hamiltonian
%   [MHz] of the spins 'Spins' of the system
%   'SpinSystem'. Spins = 1 is the first electron.
%
%   If 'sparse' is given, the matrix is returned in sparse format.
