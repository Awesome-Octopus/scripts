% nnint  Nuclear-nuclear spin interaction Hamiltonian 
%
%   Hnn = nnint(SpinSystem)
%   Hnn = nnint(SpinSystem,nucSpins)
%   Hnn = nnint(SpinSystem,nucSpins,'sparse')
%
%   Returns the nuclear-nuclear spin interaction (NNI)
%   Hamiltonian, in MHz.
%
%   Input:
%   - SpinSystem: Spin system structure. NNI
%       parameters are in the nn and nnFrame fields.
%   - nucSpins: Indices of nuclear spins for which the NNI should be
%       computed. E.g. [1 3] indicates the first and third nuclear spin.
%       If nucspins is omitted or empty, all nuclei are included.
%   - 'sparse': If given, the matrix is returned in sparse format.
%
%   Output:
%   - Hnn: Hamiltonian matrix containing the NNI for nuclear spins specified
%       in nucSpins.
