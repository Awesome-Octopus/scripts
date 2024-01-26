% crystalfield  Crystal-field Hamiltonian for the orbital angular momentu,
%
%   F = zfield(SpinSystem)
%   F = zfield(SpinSystem,OAMs)
%   F = zfield(SpinSystem,OAMs,'sparse')
%
%   Returns the crystal-field
%   Hamiltonian [MHz] of the system SpinSystem.
%
%   If the vector OAMs is given, the Crystal-field of only the
%   specified orbital angular momentums is returned (1 is the first, 2 the
%   second, etc). Otherwise, all orbital angular momentums are included.
%
%   If 'sparse' is given, the matrix is returned in sparse format.
