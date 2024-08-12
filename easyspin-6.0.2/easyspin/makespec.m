% makespec   Construct stick spectrum from peak data 
%
%  spec = makespec(Range,nPoints,Pos)
%  spec = makespec(Range,nPoints,Pos,Amp)
%  [x,spec] = makespec(...)
%
%  Constructs a stick spectrum from peak positions and intensities.
%
%   Range    limits of the abscissa [minX maxX]
%   nPoints  length of spectral vector
%   Pos      array of peak positions
%   Amp      array of peak amplitudes
%            if omitted all amplitudes are set to 1
%
%   spec     spectrum
%   x        abscissa vector
%
%   Pos and Amp must contain the same number of elements
