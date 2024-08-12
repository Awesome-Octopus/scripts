% curry  Computation of magnetometry data (magnetic moment, susceptibility)
%
%   curry(Sys,Exp)
%   curry(Sys,Exp,Opt)
%   ... = curry(...)
%
%    Calculates the magnetic moment, the molar static magnetic
%    susceptibility and related quantities of a sample for given values of
%    applied magnetic field and temperature.
%
%    Input:
%      Sys    spin system
%        .TIP          temperature-independent molar susceptibility
%                      (in SI units, m^3 mol^-1)
%
%      Exp    experimental parameter
%        .Field        list of field values (mT)
%        .Temperature  list of temperatures (K)
%
%      Opt    calculation options
%        .Output    string of keywords defining the outputs
%                   'mu'        single-center magnetic moment along field
%                   'mumol'     molar magnetic moment (magnetization) along field
%                   'muBM'      single-center magnetic moment along field,
%                                 as multiple of Bohr magnetons
%                   'mueff'     effective magnetic moment (unitless)
%                   'chi'       single-center magnetic susceptibility,
%                                 component along field
%                   'chimol'    molar magnetic susceptibility,
%                                 component along field
%                   'chimolT'   chimol times temperature
%                   '1/chimol'  inverse of chimol
%                   The default is 'muBM chimol'.
%        .Units     'SI' (for SI units, default) or 'CGS' (for CGS-emu units)
%        .Method    calculation method, 'operator' (default) or 'energies'
%        .GridSize  grid size for powder average
%        .Spins     electron spin indices, for spin-selective calculation
%        .deltaB    field step to use to calculate susceptibility (mT)
%
%
%    Output depends on the settings in Opt.Output. If Opt.Output is not given,
%    two outputs are provided:
%      muzBM    magnetic moment along zL axis, as multiple of Bohr magnetons
%                 (value is the same in SI and CGS-emu)
%      chimol   molar susceptibility, zLzL component, in SI units (m^3 mol^-1)
%
%    The size of outputs is nB x nT, where nB is the number of field values in
%    Exp.Field and nT is the number of temperature values in Exp.Temperature.
%
%   If no output argument is given, the computed data are plotted.
