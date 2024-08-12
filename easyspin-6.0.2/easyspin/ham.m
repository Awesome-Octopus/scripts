% ham  Spin Hamiltonian 
%
%   [H0,mux,muy,muz] = ham(Sys)
%   [H0,muzL] = ham(Sys,B0)
%   H = ham(Sys,B0)
%   H = ham(Sys)
%   __ = ham(__,'sparse')
%
%   Constructs the spin Hamiltonian matrix, or its field-independent and
%   field-dependent components, for the spin system Sys.
%
%   Input:
%     Sys       spin system specification structure
%     B0        vector specifying the static magnetic field (in mT) in the
%               molecular frame (e.g. [350 0 0] is along the molecular x axis)
%     'sparse'  if given, sparse instead of full matrices are returned
%
%   Output:
%     H            complete spin Hamiltonian (MHz)
%     H0           field-independent part of spin Hamiltonian (MHz)
%     mux,muy,muz  magnetic dipole moment operators (MHz/mT)
%                  along x, y, z axes of molecular frame
%     muzL         magnetic dipole moment operator (MHz/mT)
%                  along z axis of lab frame
%
%   The spin Hamiltonian components are defined such that
%           H = H0 - B0(1)*mux - B0(2)*muy - B0(3)*muz
%           H = H0 - norm(B0)*muzL
