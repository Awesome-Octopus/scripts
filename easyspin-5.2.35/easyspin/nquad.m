% nquad  Nuclear quadrupole interaction Hamiltonian
%
%   Hnq = nquad(SpinSystem)
%   Hnq = nquad(SpinSystem,Nuclei)
%   Hnq = nquad(SpinSystem,Nuclei,'sparse')
%
%   Returns the nuclear quadrupole interactions (NQI)
%   Hamiltonian, in MHz, of the spin system given in
%   SpinSystem. If the vector Nuclei is given, only
%   the NQI of the specified nuclei are computed. 1 is
%   the first nucleus, 2 the second, etc.
%
%   If 'sparse' is given, the matrix is returned in sparse format.
