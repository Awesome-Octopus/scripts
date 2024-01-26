% resonator  Resonator reflection coefficient
%
%  Gamma = resonator(nu,nu0,Qu,beta)
%
%  Calculates the reflection coefficient for a reflection resonator,
%  i.e. the fraction of voltage from an incoming wave that is
%  reflected by a resonator.
%
%  Inputs:
%    nu      array of frequency range, GHz
%    nu0     resonator frequency, GHz
%    Qu      unloaded Q-factor of the resonator
%    beta    coupling coefficient
%            beta<1: undercoupled
%            beta=1: critically coupled (matched)
%            beta>1: overcoupled
%
%    If nu is [], then an appropriate frequency range is chosen automatically.
%
%  Outputs:
%    Gamma   Voltage reflection coefficient (ratio of E-field amplitudes
%            of reflected and incoming wave). Gamma is the _voltage_ reflection
%            coefficient. The power reflection coefficient is abs(Gamma)^2.
%
%    If no output is requested, the results are plotted.
%
%  Example:
%    nu = linspace(9.2,9.8,1001);
%    nu0 = 9.5;
%    Qu = 1000;
%    beta = 1;
%    resonator(nu,nu0,Qu,beta);
