%esfit_simplex    Nelder/Mead simplex minimization algorithm
%
%   xmin = esfit_simplex(fcn,x0,lb,ub,FitOpt)
%   [xmin,info] = ...
%
%   Tries to find x that minimizes fcn(x), starting at x0. Opt are options
%   for the minimization algorithm. Any additional parameters are passed
%   to fcn, which can be a string or a function handle.
%
%   Input:
%     fcn   ... function to minimize, f(x), where is an array of parameters
%     x0    ... initial parameter values
%     lb    ... lower bounds for parameters
%     ub    ... upper bounds for parameters
%     FitOpt... structure with options
%       .delta         edge length of initial simplex
%       .SimplexPars   [rho, chi, psi, sigma]
%          rho ...     reflection coefficient
%          chi ...     expansion coefficient
%          psi ...     contraction coefficient
%          sigma .     reduction coefficient
%          The default is [1,2,1/2,1/2] for one- and two-dimensional
%          problems, and adaptive for higher dimensions.
%       .maxTime       maximum time allowed, in minutes
%
%   Output:
%     xmin  ... parameter vector with values of best fit
%     info  ... structure with additional information (initial simplex,
%               last simplex, number of iterations, time elapsed)
