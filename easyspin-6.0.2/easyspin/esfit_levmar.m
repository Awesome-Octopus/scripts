%esfit_levmar  Levenberg-Marquardt nonlinear least squares fitting
%
%   xfit = esfit_levmar(fcn,x0,lb,ub)
%   ... = esfit_levmar(fcn,x0,lb,ub,FitOpt)
%   [xfit,info] = ...
%
%   Tries to find x that minimizes sum(fcn(x).^2), starting at x0.
%   fcn is a function that provides a vector of residuals.
%
% Input
%   fcn    residuals vector
%   x0     starting vector in parameter space
%   lb     lower bounds of parameters
%   ub     upper bounds of parameters
%   FitOpt structure with algorithm parameters
%     .lambda    starting value of Marquardt parameter
%     .Gradient  termination threshold for gradient
%     .TolStep   termination threshold for step
%     .maxTime   termination threshold for time
%     .delta     step width for Jacobian approximation
%     .Verbosity print detail level
%     .IterFcn   function that is called after each iteration
%
% Output
%   xfit    Converged vector in parameter space
%   info    structure with fitting information
%            .F  function value at minimum
%            .norm_g
%            .norm_h
%            .Je Jacobian  estimate
%            .lambda
%            .nIterations  number of interations
%            .stop
%            .nEvals       number of function evaluations
%            .msg          message
