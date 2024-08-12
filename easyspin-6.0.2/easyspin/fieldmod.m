% fieldmod  field modulation 
%
%   yMod = fieldmod(B,spc,ModAmp);
%   yMod = fieldmod(B,spc,ModAmp,Harmonic);
%   fieldmod(...)
%
%   Computes the effect of field modulation on an EPR absorption spectrum.
%
%   Input:
%   - B: magnetic field axis vector, mT
%   - spc: absorption spectrum
%   - ModAmp: peak-to-peak modulation amplitude, mT
%   - Harmonic: harmonic (0, 1, 2, ...); default is 1
%
%   Output:
%   - spcMod: pseudo-modulated spectrum
%
%   If no output variable is given, fieldmod plots the
%   original and the modulated spectrum.
%
%   Example:
%
%     B = linspace(300,400,1001);  % mT
%     spc = lorentzian(B,342,4);
%     fieldmod(B,spc,20);
