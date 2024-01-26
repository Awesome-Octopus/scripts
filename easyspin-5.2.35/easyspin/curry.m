% curry  Computation of magnetometry data
%         (magnetic moment, susceptibility)
%
%   curry(Sys,Exp)
%   curry(Sys,Exp,Opt)
%   muz = curry(...)
%   [muz,chizz] = curry(...)
%
%    Calculates the magnetic moment and the molar static magnetic
%    susceptibility of a powder sample for given values of
%    applied magnetic field and temperature.
%
%    Input:
%      Sys    spin system
%      Exp    experimental parameter
%        Field        list of field values (mT)
%        Temperature  list of temperatures (K)
%        
%        mField, mTemperature 
%        chiField, chiTemperature
%                     list of field values and temperatures for the
%                     calculation of the magnetic moment (mField &
%                     mTemperature) and the magnetic susceptibility
%                     (chiField & chiTemperature)
%
%
%      Opt    calculation options
%        nKnots       number of knots for powder average
%        Output       list of keywords defining the output and the order of
%                     it
%                     the following keywords are allowed:
%                     'MvsB', 'MvsBCGS', 'MvsBSI', 'Chi',
%                     'ChiT', '1overChi', 'MuEff','ChiSI'
%                     'ChiTSI', '1overChiSI', 'MuEffSI',
%                     'ChiCGS', 'ChiTCGS', '1overChiCGS', 
%                     'MuEffCGS', 'OneColoumn'
%        Method       'operator' (default) or 'energies'
%                     calculation method
%  
%                     
%    Output (if Opt.Output is not given):
%      muz      magnetic moment along zL axis
%               in units of Bohr magnetons
%      chizz    molar susceptibility, zLzL component
%               in SI units (m^3 mol^-1)
%    
%    zL is the direction of the applied static magnetic field
%
%    The size of muz and chizz is nB x nT, where nB is the number of
%    field values in Exp.Field and nT is the number of temperature values
%    in Exp.Temperature.
%
%   If no output argument is given, the computed data are plotted.
