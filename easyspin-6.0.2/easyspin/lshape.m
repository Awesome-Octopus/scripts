% lshape  Gaussian and Lorentzian line shapes 
%
%   y = lshape(x,x0,fwhm)
%   y = lshape(x,x0,fwhm,diff)
%   y = lshape(x,x0,fwhm,diff,alpha)
%   y = lshape(x,x0,fwhm,diff,alpha,phase)
%
%   General normalized line shape function, a linear
%   combination of Lorentzian and Gaussian lineshapes.
%
%   Input:
%   - x: Abscissa vector
%   - x0: Center of the linshape function
%   - fwhm: Full width at half height. If alpha is
%     given, and Gaussian and Lorentzian components should
%     have different fwhm, specify a vector
%         [fwhmGauss fwhmLorentz]
%   - diff: Derivative. 0 is no derivative, 1 first,
%     2 second, -1 the integral with -infinity as lower
%     limit. If omitted, 0 is the default.
%   - alpha: Weighting factor for linear combination
%       alpha*Gaussian + (1-alpha)*Lorentzian.
%     Gaussian is default.
%   - phase: phase (0 pure absorption, pi/2 pure dispersion)
%
%   Output:
%   - y: Vector of function values for abscissa x
