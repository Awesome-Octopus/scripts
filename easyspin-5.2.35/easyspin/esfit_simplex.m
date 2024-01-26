%esfit_simplex    Nelder/Mead downhill simplex minimization algorithm
%
%   xmin = esfit_simplex(fcn,x0,FitOpt,...)
%   [xmin,info] = ...
%
%   Tries to minimize fcn(x), starting at x0. FitOpt are options for
%   the minimization algorithm. Any additional parameters are passed
%   to fcn, which can be a string or a function handle.
%
%   Options:
%     FitOpt.SimplexPars   [rho chi psi sigma]
%                          rho ... reflection coefficient, default 1
%                          chi ... expansion coefficient, default 2
%                          psi ... contraction coefficient, default 0.5
%                          sigma . reduction coefficient, defautlt 0.5
%     FitOpt.maxTime       maximum time allowed, in minutes
%     FitOpt.delta         edge length of initial simplex
%
%   Output:
%     xmin  ... parameter vector with values of best fit
%     info  ... structure with additional information (initial simplex,
%               last simplex, number of iterations, time elapsed)
