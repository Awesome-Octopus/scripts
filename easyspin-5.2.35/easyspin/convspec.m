% convspec  Spectrum convolution with line shapes 
%
%   out = convspec(spec,df,fwhm)
%   out = convspec(spec,df,fwhm,deriv)
%   out = convspec(spec,df,fwhm,deriv,alpha)
%   out = convspec(spec,df,fwhm,deriv,alpha,phase)
%
%   Convolutes the spectral data array spec with a line
%   shape. df specifies the abscissa step, fwhm the
%   fwhm line width, deriv is the derivative (default 0),
%   alpha the shape parameter (1 is Gaussian, 0 is Loren-
%   tzian, default 1, between 0 and 1 a pseudovoigtian is
%   computed). phase is the phase of the Lorentzian component
%   (0 pure absorption, pi/2 pure dispersion).
%   out is the convoluted spectrum.
%
%   If the spectrum is more than 1D, different parameters
%   for each dimension can be defined. For 2D, e.g. df = 2
%   means that the abscissa step is 2 for both dimensions.
%   df = [2 3] means it is 2 for dim 1 and 3 for dim 2.
%   The parameters fwhm, deriv and alpha work similarly.
%   If fwhm is 0 for a certain dimension, convolution
%   along this dimension will be skipped.
%
%   Examples
%     spec = zeros(101,101);
%     spec(51,51) = 1; spec(30,60) = 3;
%     w = convspec(spec,.1,1);
%     pcolor(w);
%
%    convolutes in both dimensions with a Gaussian
%    with fwhm=1. The abscissa step is 0.1 for both
%    dimensions.
%
%     spec = zeros(100,100);
%     spec(51,51) = 1; spec(30,60) = 3;
%     w = convspec(spec,1,[3 6],[0 1],[0 1]);
%     surf(w);
%
%    convolutes in dim 1 with a Lorentzian fwhm=5 and
%    in dim 2 with a first derivative Gaussian fwhm=3.
