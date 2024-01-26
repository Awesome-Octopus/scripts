% mhz2mt   Conversion from MHz to mT
%
%   x_mT = mhz2mt(x_MHz)
%   x_mT = mhz2mt(x_MHz,g)
%
%   Converts the value in x_MHz, assumed to be in units
%   of MHz (Megahertz) to mT (Millitesla), returning
%   the result in x_mT.
%
%   For the conversion, the g factor given as second
%   parameter is used. If it is not given, the g factor
%   of the free electron (gfree) is used.
%
%   x_MHz can be a vector of values. In this case, g
%   can be a scalar or a vector of the same size as x_MHz.
