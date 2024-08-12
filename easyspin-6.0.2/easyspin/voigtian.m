% voigtian  Voigtian line shape function 
%
%   ya = voigtian(x,x0,fwhmGL)
%   ya = voigtian(x,x0,fwhmGL,diff)
%   ya = voigtian(x,x0,fwhmGL,diff,phase)
%   [ya,yd] = voigtian(...)
%
%    Computes an area-normalized Voigt line shape function, which is
%    the convolution of a Gaussian and a Lorentzian line shape function.
%
%   Input:
%   - x:     abscissa vector
%   - x0:    center of the lineshape function
%   - fwhmGL:  full widths at half height [fwhm_Gauss fwhm_Lorentz]
%   - diff:  derivative. 0 is no derivative, 1 first,
%            2 second and so on, -1 the integral with -infinity
%            as lower limit. 0 is the default.
%   - phase: phase rotation, mixes absorption and dispersion.
%            phase=pi/2 puts dispersion signal into ya
%
%   Output:
%   - ya:   absorption function values for abscissa x
%   - yd:   dispersion function values for abscissa x
