% ham_cf  Crystal-field Hamiltonian for the orbital angular momentum
%
%   H = ham_cf(Sys)
%   H = ham_cf(Sys,OAMs)
%   H = ham_cf(Sys,OAMs,'sparse')
%
%   Returns the crystal-field Hamiltonian (in MHz) of the spin system
%   Sys, utilizing the fields Sys.L, Sys.CF1, Sys.CF2, etc.
%
%   If the vector OAMs is given, the crystal field of only the
%   specified orbital angular momenta is returned (1 is the first in Sys.L,
%   2 the second, etc). Otherwise, all orbital angular momenta are included.
%
%   If 'sparse' is given, the matrix is returned in sparse format.
