% lorentzian  Lorentzian line shape 
%
%   ya = lorentzian(x,x0,fwhm)
%   ya = lorentzian(x,x0,fwhm,diff)
%   ya = lorentzian(x,x0,fwhm,diff,phase)
%   [ya,yd] = lorentzian(...)
%
%   Computes area-normalized Lorentzian absorption and dispersion
%   line shapes or their derivatives.
%
%   Input:
%   - x:     abscissa vector
%   - x0:    center of the lineshape function
%   - fwhm:  full width at half height
%   - diff:  derivative. 0 is no derivative, 1 first,
%            2 second and so on, -1 the integral with -infinity
%            as lower limit. 0 is the default.
%   - phase: phase rotation, mixes absorption and dispersion.
%            phase=pi/2 puts dispersion signal into ya
%
%   Output:
%   - ya:    absorption function values for abscissa x
%   - yd:    dispersion function values for abscissa x
