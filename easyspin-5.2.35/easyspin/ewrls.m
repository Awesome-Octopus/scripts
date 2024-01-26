% ewrls   Exponentially weighted recursive least squares adaptive filter
%
%   y = ewrls(data,p,lambda)
%   y = ewrls(data,p,lambda,nPreAvg)
%   y = ewrls(data,p,lambda,nPreAvg,delta)
%   y = ewrls(data,p,lambda,nPreAvg,delta,dir)
%   ewrls(...)
%
%   Given a series of scans in the 2D data array data, computes a filtered
%   average using the exponentially weighted recursive least squares (RLS)
%   adaptive filter method.
%   If you omit the output (y), ewrls plots the result.
%
%     data     2D data array, each column representing one scan
%     p        filter length, typically 5-50
%     lambda   memory factor, typically 0.96-0.99
%     nPreAvg  number of scans to obtain desired signal
%              if 0, all scans are averaged. default 0
%     delta    regularization parameter (default 100)
%     dir      filtering direction (default 'fb')
%              'f' forward, 'b' backward, 'fb' average of both
%     y        filtered averaged scans, a single column
%
%   Example:
%     y = ewrls(data,32,0.96)
