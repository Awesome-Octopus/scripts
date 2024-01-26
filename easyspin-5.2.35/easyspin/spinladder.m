% spinladder    Manifold Hamiltonian of a spin ladder
%
%    CSys = spinladder(Sys)
%    [CSys,En] = spinladder(Sys)
%    ... = spinladder(Sys,Temp)
%    spinladder(...)
%
%    Given a exchange-coupled two-electron-spin system in Sys,
%    this function computes the spin Hamiltonians for the
%    various spin manifolds in the coupled representation,
%    incl. g, A and D tensors, assuming the strong-exchange
%    limit.
%
%    CSys is a cell array that contains the coupled-spin
%         systems sorted by energy.
%    En   contains the center-of-gravity energies as
%         determined by the exchange coupling.
%
%    If no output is requested, spinladder prints some
%    information about the coupled manifolds.
%
%    If a temperature (in kelvin) is given, populations for
%    each manifold are computed and returned in the fields
%    CSys{}.weight.
