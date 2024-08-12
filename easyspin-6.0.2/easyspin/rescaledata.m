% rescaledata    Rescaling of data
%
%   yscaled = rescaledata(y,mode)
%   yscaled = rescaledata(y,yref,mode)
%   yscaled = rescaledata(...,region)
%   [yscaled, scale] = rescaledata(...)
%
% Rescales the data in array y. If given, yref serves
% as the reference. mode determines how the scale factor is calculated.
% Positive scaling is enforced, i.e. the rescaled data is never inverted.
%
% Inputs:
%   y           data to be rescaled, 1D array
%   yref        reference data (used by 'maxabs' and 'lsq'), 1D array
%   mode:
%     'maxabs'  scales y such that maximum magnitude of yscaled is 1 (if
%               yref is not given) or max(abs(yref)) (if yref is given)
%     'lsq'     least-squares fit of a*y to yref; yref is needed
%     'int'     normalize integral (sum of datapoints) to 1
%     'dint'    normalize double integral (sum of cumsum of datapoints) to 1
%     'none'    no scaling
%
% Outputs:
%   yscaled   rescaled data vector
%   scale     scaling factor
