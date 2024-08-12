% ham_hf  Hyperfine interaction Hamiltonian 
%
%   Hhf = ham_hf(System)
%   Hhf = ham_hf(System,eSpins)
%   Hhf = ham_hf(System,eSpins,nSpins)
%   Hhf = ham_hf(System,eSpins,nSpins,'sparse')
%
%   Returns the hyperfine interaction Hamiltonian (in units of MHz) between
%   electron spins 'eSpins' and nuclear spins 'nSpins' of the system
%   'System'. eSpins=1 is the first electron spins, nSpins=1 is the first
%   nuclear spin.
%
%   If 'sparse' is given, the matrix is returned in sparse format.
