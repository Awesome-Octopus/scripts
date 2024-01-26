% hilberttrans  Hilbert transform
%
%   yh = hilberttrans(y)
%
%   Computes the Hilbert transform of input vector y. It can
%   be used to compute a dispersion signal from an absorption
%   signal.
%
%   Input:
%     - y:    input vector (spectrum)
%   Output:
%     - yh:   Hilbert transform
%               real part: identical to y
%               imaginary part: actual Hilbert transform of y
%
%   Example:
%      x = linspace(-1,1,1000);
%      yabs = gaussian(x,0,0.1,0);
%      ydisp = hilberttrans(yabs);
%      ydisp = imag(ydisp);
%      plot(x,yabs,x,ydisp)
