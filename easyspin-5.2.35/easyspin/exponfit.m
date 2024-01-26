% exponfit  Fits exponential(s) to 1D and 2D data
%
%   [k,c] = exponfit(x,y)
%   ...   = exponfit(x,y,nExps)
%   ...   = exponfit(x,y,nExps,'noconst')
%   [k,c,yFit] = exponfit(...)
%   k = exponfit(...)
%
%   Computes least square fits for single or double exponential
%   decays or recoveries. For matrices it works along columns and
%   returns the result in corresponding matrix form. The fitting model
%   is of the form
%
%      c(1) + c(2)*exp(-k(1)*x) + c(3)*exp(-k(2)*x)
%
%   Input:
%     x          abscissa values, a vector
%     y          ordinate values, a vector or a matrix
%     nExps      number of exponents to include (1 or 2, default is 1)
%     'noconst'  if given, c(1) is set to 0 and not included in the fit
%
%   Output:
%     k       exponents
%     c       linear coefficients
%     yFit    the fitting model function
